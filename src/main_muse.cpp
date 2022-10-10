/* 
Statement of MuSE 2.0
MuSE 2.0, powered by parallel computing by taking advantage of multi-core resource on a machine, 
and a more efficient way of memory allocation, is about 40-50x faster than MuSE 1.0 on mutation 
calling with whole exome sequencing or whole genome sequencing data. We thank Mehrzad Samadi and 
his team from Nvidia Corporation, including Tong Zhu, Timothy Harkins and Ankit Sethia, for their 
contributions of implementing accelerating techniques in the ‘MuSE call’ step.
*/

#include <iostream>
#include <cstring>
#include <sys/stat.h>
#include "htslib/bgzf.h"
#include "muse_reader.h"
#include "timer.h"
#include "tabix.h"
#include <omp.h>

using namespace std;

int tid_global = 0;
int64_t pos_global = -1;

void muse_sump(const char *inFile, const char *outFile, const char *dbsnpFile, bool isWGS, bool isWES, int num_threads, int argc, char *argv[]);

void monitorFun(PBReader* reader, std::atomic<uint32_t>& processQSize, PileupSpscQ& writeQ, std::atomic<bool>& monitorFlag){
	while(monitorFlag.load()){
		std::time_t result = std::time(nullptr);
		std::string timeNowStr(std::asctime(std::localtime(&result)));
		std::string timeNow(timeNowStr, 11, 8);
		fprintf(stderr, "[%s]\t%s:%ld\n", timeNow.c_str(), reader->chrName(tid_global), pos_global);
		fprintf(stderr, "BamRead %d processQSize %d writeQSize %lu readPool %lu\n", reader->getReadQueueSize(), processQSize.load(), writeQ.read_available(), reader->getReadPoolSize());
		std::this_thread::sleep_for(std::chrono::seconds(1));
	}
}

void readerFun(PBReader* reader){
	reader->fillBamQueue();		
	reader->setDone();
}

uint32_t bam_calend(bam1_t* b)
{	auto& c = b->core;
	auto cigar = bam_get_cigar(b);
	uint32_t k, end;
	end = c.pos;
	for (k = 0; k < c.n_cigar; ++k) {
		int op = cigar[k] & BAM_CIGAR_MASK;
		if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP)
			end += cigar[k] >> BAM_CIGAR_SHIFT;
	}
	return end;
}

void processReadFun(PBReader* reader, std::atomic<bool>& gatherDone, PileUpLockFreeQ& processQ, std::atomic<uint32_t>& processQSize, PileupSpscQ& writeQ, mplp_conf_t* conf){

	bam_plp iter(conf->max_depth);
	bam1_t* read = nullptr;
	Region* region;
	while(read = reader->getNextFilteredRecordFromQueue()){
		node Node = {read, bam_calend(read)};
		region = iter.push(Node);
		if (region){
			while(processQSize.load() > PILEUP_SIZE);
			assert(processQ.push(region));
			++processQSize;
			while(writeQ.push(region) == false);
		}
	}

	while(region = iter.pop()){
		while(processQSize.load() > PILEUP_SIZE);
		assert(processQ.push(region));
		++processQSize;
		while(writeQ.push(region) == false);
	}
	gatherDone.store(true);

}

void pushToReadPool(Region* localPileup, PBReader* reader){
	for (auto& x: localPileup->nodes){
		if (x.end - 1 <= localPileup->end){
			reader->release(x.b);
			x.b = nullptr;
		}
	}
	localPileup->nodes.clear();
}

void writeFun(std::atomic<bool>& gatherDone, PileupSpscQ& writeQ, string& outputName, PBReader* reader, mplp_conf_t* conf){
	outputName += ".MuSE.txt";
	PBLocalFile* file = new PBLocalFile(outputName.c_str(), "w");
	if (file->isFailed()){
		cerr << "Cannot write to output file. Exiting..." << endl;
		exit(EXIT_FAILURE);
	}
	WriteHeader(file, conf->argc, conf->argv, reader->mHeader_tumor, reader->mHeader_normal, conf->ref->fileName);

	Region* localPileup;

	while(!gatherDone.load() || writeQ.read_available()){
		if (!writeQ.pop(localPileup))
			continue;
		tid_global = localPileup->tid;
		pos_global = localPileup->start;
		while(localPileup->ifDone.load() == false);

		file->pbwrite(localPileup->output.c_str(), 1, localPileup->output.size());

		pushToReadPool(localPileup, reader);
		delete localPileup;
	}

	delete file;
}

void mpileup(mplp_conf_t& conf, string& tumorName, string& normalName, string& outputName, int threadNum){
	construct_convert_table();
	PileUpLockFreeQ processQ(PILEUP_SIZE);
	std::atomic<uint32_t> processQSize(0);
	PileupSpscQ writeQ;

	std::atomic<bool> gatherDone(false);
	std::atomic<bool> monitorFlag(true);
	PBReader reader(tumorName, normalName, &conf);
	PBThread moniterThread("monitor", monitorFun, &reader, std::ref(processQSize), std::ref(writeQ), std::ref(monitorFlag));
	PBThread readerThread("reader", readerFun, &reader);

	PBThread readProcessThread("processReads", processReadFun, &reader, std::ref(gatherDone), std::ref(processQ), std::ref(processQSize), std::ref(writeQ), &conf);
	
	vector<PBThread> processThreads(threadNum);
	for (int i = 0; i < threadNum; ++i){
		processThreads[i] = PBThread("worker_" + to_string(i), processPileup, &reader,  &conf, std::ref(gatherDone), std::ref(processQ), std::ref(processQSize));
	}
	PBThread writerThread("wirter", writeFun, std::ref(gatherDone), std::ref(writeQ), std::ref(outputName), &reader, &conf);
	readerThread.join();

	readProcessThread.join();
	for (auto& x: processThreads)
		x.join();
	writerThread.join();
	monitorFlag.store(false);

	moniterThread.join();
}

void get_MuseCallOpts(int argc, char* argv[]){
    
    int c;
    string tumorName, normalName, outputName, refName, threadNum_c = "1";
    int threadNum = 1;

	while((c = getopt(argc, argv, "f:n:O:")) >= 0) {
		switch(c) {
			case 'f': refName           = optarg;         break;
			case 'n': threadNum_c       = optarg;         break;
			case 'O': outputName        = optarg;         break;
		}
	}

    if (argc == 2){
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   MuSE call [options] tumor.bam matched_normal.bam\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "         -f FILE    faidx indexed reference sequence file\n");
		fprintf(stderr, "         -O STR     output file name (suffix '.MuSE.txt' is\n");
		fprintf(stderr, "                    automatically added)\n");
		fprintf(stderr, "         -n INT     number of cores specified (default=1)\n");
		fprintf(stderr, "\n");

        exit(EXIT_FAILURE);
    }

    if (argc <= optind + 2){
        cerr << "Wrong argument. Existing..." << endl;
        exit(EXIT_FAILURE);
    }
    
    try{
        threadNum = stoi (threadNum_c);
    }
    catch(const std::exception& e){
        cerr << e.what() << endl;
        exit(EXIT_FAILURE);
    }

    if (threadNum < 1){
        cerr << "Number of cores cannot be less than 1. Exiting..." << endl;
		exit(EXIT_FAILURE);
    }

	if (refName.empty()){
		cerr << "No reference file specified. Exiting..." << endl;
		exit(EXIT_FAILURE);
	}

	if (outputName.empty()){
		cerr << "No output file name specified. Exiting..." << endl;
		exit(EXIT_FAILURE);
	}

    tumorName = argv[argc-2];
    normalName = argv[argc-1];

    mplp_conf_t mplp(refName, argc, argv);
    mpileup(mplp, tumorName, normalName, outputName, threadNum);

}

void get_MuseSumpOpts(int argc, char *argv[]){
    int        c;
    const char *outFile     = NULL;
    const char *inFile      = NULL;
    const char *dbsnpFile   = NULL;
    bool       isWGS        = false;
    bool       isWES        = false;

    const char *threadNum_c = "0";
    int threadNum;
    // command options
    //

	while((c = getopt(argc, argv, "I:O:D:n:GE")) >= 0) {
        switch(c) {
        case 'I': inFile    = optarg; break;
        case 'O': outFile   = optarg; break;
        case 'D': dbsnpFile = optarg; break;
	case 'n': threadNum_c  = optarg; break; 
        case 'G': isWGS     = true;   break;
        case 'E': isWES     = true;   break;
        }
    }
    if(argc == 1) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   MuSE sump [options]\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "         -I FILE    single input file generated by 'MuSE call'\n");
        fprintf(stderr, "         -G         input generated from whole genome sequencing data\n");
        fprintf(stderr, "         -E         input generated from whole exome sequencing data\n");
        fprintf(stderr, "         -O STR     output file name (VCF format)\n");
	fprintf(stderr, "         -n int     number of cores specified (default=1)\n");
        fprintf(stderr, "         -D FILE    dbSNP vcf file that should be bgzip compressed,\n");
        fprintf(stderr, "                    tabix indexed and based on the same reference\n");
        fprintf(stderr, "                    genome used in 'MuSE call'\n");
        fprintf(stderr, "\n");
        exit(EXIT_FAILURE);
    }

    // must identify the data type but cannot be both
    //
    if(isWGS && isWES) {
        cerr << "Option -G and -E cannot be selected at the same time. Exiting..." << endl;
        exit(EXIT_FAILURE);
    }
    if((!isWGS) && (!isWES)) {
        cerr << "Must identify the sequencing data type using either option -G or option -E. Exiting..." << endl;
        exit(EXIT_FAILURE);
    }

    // check input and output files
    //
    if(!inFile) {
        cerr << "Please designate the input file. Exiting..." << endl;
        exit(EXIT_FAILURE);
    }
    if(!outFile) {
        std::cerr << "Please designate the output file name. Exiting..." << endl;
        exit(EXIT_FAILURE);
    }

    try{
	threadNum = stoi (threadNum_c);
    }
    catch(const std::exception& e){
	cerr << e.what() << endl;
	exit(EXIT_FAILURE);
    }

    if (threadNum < 1){
	cerr << "Number of cores cannot be less than 1. Exiting..." << endl;
	exit(EXIT_FAILURE);
    } 
    
    // check if dbSNP file was bgzipped
    //
	if(dbsnpFile) {
        BGZF *bgzf_handle = bgzf_open(dbsnpFile, "r");
        if (bgzf_handle){
            if(bgzf_compression(bgzf_handle) != 2){
                fprintf(stderr,"%s cannot be open.\n", dbsnpFile);
                bgzf_close(bgzf_handle);
                exit(EXIT_FAILURE);
            }
        }else{
            fprintf(stderr,"Was bgzip used to compress %s?\n", dbsnpFile);
            bgzf_close(bgzf_handle);
            exit(EXIT_FAILURE);
        }
    }

	// check the index of dbSNP file
    // Common source of errors: new VCF is used with an old index.  On some systems, stat on non-existent files returns undefined value for sm_mtime
    //
    if(dbsnpFile) {
        struct stat stat_tbi,stat_vcf;
        char *fnidx = (char *)calloc(strlen(dbsnpFile) + 5, 1);
        strcat(strcpy(fnidx, dbsnpFile), ".tbi");
        stat(fnidx, &stat_tbi);
        stat(dbsnpFile, &stat_vcf);
        if(stat_vcf.st_mtime > stat_tbi.st_mtime) {
            fprintf(stderr, "The index of %s either does not exist or is older than the vcf file.\n", dbsnpFile);
            free(fnidx);
            exit(EXIT_FAILURE);
        }
        free(fnidx);
    }

    int num_threads = 1;

#ifdef _OPENMP
#pragma omp parallel
    {
#pragma omp master
        num_threads = omp_get_num_threads();
    }
#else
#endif

    if (threadNum > num_threads) threadNum = num_threads;
    
#ifdef _OPENMP 
    omp_set_num_threads(threadNum);
#endif
    
    muse_sump(inFile, outFile, dbsnpFile, isWGS, isWES, threadNum, argc, argv);
}

//================================================================================================= Main
static int usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: MuSE (Tools for calling somatic point mutations)\n\n");
	fprintf(stderr, "Version: %s\n", Version.c_str());
    fprintf(stderr, "         Build Date %s\n", buildDate.c_str());
    fprintf(stderr, "         Build Time %s\n\n", buildTime.c_str());
    fprintf(stderr, "Usage:   MuSE <command> [options]\n\n");
    fprintf(stderr, "Command: call    Call somatic point mutations.\n");
    fprintf(stderr, "                 Raw output file will be created.\n");
    fprintf(stderr, "         sump    Generate final calls in VCF format.\n");
    fprintf(stderr, "                 Besides PASS calls, tier-based calls are\n");
    fprintf(stderr, "                 reported using dynamic cutoff searching\n");
    fprintf(stderr, "                 approach. If there are multiple raw output\n");
    fprintf(stderr, "                 files from WGS data, please have them\n");
    fprintf(stderr, "                 concatenated first.\n");
    fprintf(stderr, "\n");
    exit(EXIT_FAILURE);
}


int main(int argc, char* argv[]){

    if(argc < 2)
        return usage();

	if(strcmp(argv[1], "call") == 0){
        get_MuseCallOpts(argc, argv);
	}
	else if(strcmp(argv[1], "sump") == 0) {
		get_MuseSumpOpts(argc-1, argv+1);
	}
    else {
        fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
        return 1;
    }
	return 0;
}
