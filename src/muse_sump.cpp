#include <iostream>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <cstring>
#include <sys/stat.h>
#include <math.h>
#include <float.h>
#include <omp.h>
#include "htslib/bgzf.h"
#include "htslib/faidx.h"
#include "muse_reader.h"
#include "timer.h"
#include "tabix.h"
#include "statistics.h"
#include <chrono>

using namespace std;


struct callInfo {
	std::string ref;
	std::string alt;
	double tAltF;
	int readGroupCount;
	bool isDBSNP;
	std::string rsID;
	bool useSecondAllele;
	int tumorRefCount;
	int tumorAltCount;
	int normalRefCount;
	int normalAltCount;
	std::string altString;
	std::string tumorFormat;
	std::string normalFormat;
};

// for beta fit
//
std::vector<double> altf;
std::vector<double> betaShape(2, 2.0);
double min_diff_beta_lk = 1.E-04;
double beta_lnlike;

// for GMM fit
//
double finalMu[2];
double finalSigma[2];
double finalLambda[2];
double gmmCutoff;
double min_diff_misclassification_prob = 1.E-04;
double negative_misclassification_prob;
struct emEst {
	double loglik;
	double mu1;
	double mu2;
	double sigma1;
	double sigma2;
	double lambda1;
	double lambda2;
};
std::vector<emEst> EMEstimate;


// for data parsing 

struct dataEntry {
    int _index;
    std::string _Chrm;
    std::string _Position;
    std::string _Ref;
    std::string _Alt;
    double _TaltF;
    int _ReadGroupCount;
    bool _UseSecondAllele;
    int _TumorRefCount;
    int _TumorAltCount;
    int _NormalRefCount;
    int _NormalAltCount;
    std::string _AltString;
    std::string _TumorFormat;
    std::string _NormalFormat;
};


//================================================================================================= Utilities
std::string IntToString(int x) {
	std::ostringstream o;
    o << x;
    return o.str();
}

void muse_sump(const char *inFile, 
               const char *outFile, 
               const char *dbsnpFile, 
               bool isWGS, 
               bool isWES, 
               int num_threads, 
               int argc, 
               char *argv[]){
    
    int        minDepth     = 8;
    double     passVAF      = 0.02;
    double     baseVAF      = 0.01;
    double     vafRatio     = 0.05;
    double     minTumorVAF  = 0.005;

    int i, j, ithReplicate;
    

    // start ...
    //

    
    std::ifstream inf;
    inf.open(inFile);
    if(!inf) {
        std::cerr << "Failed to open " << inFile << '\n';
        exit(EXIT_FAILURE);
    }
    
    
    // load the input file content and remove unused lines
    //

    std::vector<std::string> Header;
    
    std::string buffer;

    std::ifstream in_f(inFile);
    in_f.seekg(0, std::ios::end);
    buffer.resize(in_f.tellg());
    in_f.seekg(0);
    in_f.read((char*)buffer.data(), buffer.size());
    
    std::vector<size_t> newline_pos;

    size_t pos = buffer.find('\n');
    while( pos != std::string::npos){
	newline_pos.push_back(pos);
	pos = buffer.find('\n', pos + 1);
    }
    
    std::string line = "";
    int snv_start_line = 0;
    for(size_t index = 0; index < newline_pos.size(); index++){
	line = buffer.substr(newline_pos[index], newline_pos[index + 1] - newline_pos[index]);
	size_t startpos = line.find_first_not_of("\n");
	if( string::npos != startpos)
	    line = line.substr(startpos);
	else
	    line = "";
	
	if(line != ""){
	    if (line.substr(0, 1) == "#") {
		Header.push_back(line);
		snv_start_line += 1;
	    }
	    else{
		break;
	    }
	}
    }

    if(snv_start_line >= newline_pos.size()){
	cerr << "MuSE.txt does not have variant entries. Exiting..." << endl;
	exit(EXIT_FAILURE);
    }
    
    //newline_pos.erase(newline_pos.begin(), newline_pos.begin() + snv_start_line);    

    int original_data_size = (int)newline_pos.size() - snv_start_line;
    
    std::vector<pair<int, int> > random_data_indices; 
    random_data_indices.reserve(original_data_size);

    std::vector<dataEntry> dataParser;
    dataParser.reserve(original_data_size);
    
#pragma omp parallel for private (line)
    for(int index = snv_start_line; index < snv_start_line + original_data_size; index++){
	line = buffer.substr(newline_pos[index], newline_pos[index + 1] - newline_pos[index]);
	size_t startpos = line.find_first_not_of("\n");
	if( string::npos != startpos)
	    line = line.substr(startpos);
	else
	    line = "";

	if(line != ""){

	    std::istringstream iss(line); 
	    std::vector<std::string> elements((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());

	    if(elements.size() != INPUT_COLUMN_COUNT) {
		std::cerr << "Incomplete line of the input file.\n";
		exit(EXIT_FAILURE);
	    }

	    int tDP = atoi(elements[4].c_str());
	    int nDP = atoi(elements[5].c_str());

	    if(tDP>=minDepth && nDP>=minDepth) {
		char ref = elements[2][0];
		char alt = elements[3][0];
		bool nAltfIsNA = true;
		double nAltF = 0.0;
		if(elements[14]!="NA") {
		    nAltfIsNA = false;
		    switch(alt) {
		    case 'A':
			nAltF = 2.0*atof(elements[14].c_str());
			break;
		    case 'C':
			nAltF = 2.0*atof(elements[15].c_str());
			break;
		    case 'G':
			nAltF = 2.0*atof(elements[16].c_str());
			break;
		    case 'T':
			nAltF = 2.0*atof(elements[17].c_str());
			break;
		    }
		}

		double tRefF = 0.0;
		switch(ref) {
		case 'A':
		    tRefF = atof(elements[7].c_str());
		    break;
		case 'C':
		    tRefF = atof(elements[8].c_str());
		    break;
		case 'G':
		    tRefF = atof(elements[9].c_str());
		    break;
		case 'T':
		    tRefF = atof(elements[10].c_str());
		    break;
		}

		double tAltF = 0.0;
		switch(alt) {
		case 'A':
		    tAltF = 2.0*atof(elements[7].c_str());
		    break;
		case 'C':
		    tAltF = 2.0*atof(elements[8].c_str());
		    break;
		case 'G':
		    tAltF = 2.0*atof(elements[9].c_str());
		    break;
		case 'T':
		    tAltF = 2.0*atof(elements[10].c_str());
		    break;
		}

		if(nAltF/tAltF <= vafRatio) {
		    bool tAltfIsDorminant = false;
		    if(0.5*tAltF > 1.0 - tRefF - 0.5*tAltF)
			tAltfIsDorminant = true;
		    if(tAltfIsDorminant) {
			
			int nRG       = atoi(elements[21].c_str());
			int tRefCount = atoi(elements[23].c_str());
			int tAltCount = atoi(elements[24].c_str());
			int nRefCount = atoi(elements[25].c_str());
			int nAltCount = atoi(elements[26].c_str());

			dataEntry _dataEntry;
			_dataEntry._index = index;
			_dataEntry._Chrm = elements[0];
			_dataEntry._Position = elements[1];;
			_dataEntry._Ref = elements[2];
			_dataEntry._Alt = elements[3] ;
			_dataEntry._TaltF = tAltF;
			_dataEntry._ReadGroupCount = nRG;

			if(elements[22]=="Y")
			    _dataEntry._UseSecondAllele = true;
			else
			    _dataEntry._UseSecondAllele = false;
			_dataEntry._TumorRefCount = tRefCount;
			_dataEntry._TumorAltCount = tAltCount;
			_dataEntry._NormalRefCount = nRefCount;
			_dataEntry._NormalAltCount = nAltCount;
			_dataEntry._AltString = elements[27];
			_dataEntry._TumorFormat = elements[28];
			_dataEntry._NormalFormat = elements[29]; 
			
#pragma omp critical
			{
			    dataParser.push_back(_dataEntry);
			}
		    }
		}
	    }
	}
    }

    buffer = "";
    int dataSize = dataParser.size();
    
    for(int index = 0; index < dataSize; index++){
	random_data_indices.push_back(make_pair(dataParser[index]._index, index));
    }

    sort(random_data_indices.begin(), random_data_indices.end());

    std::vector<std::string> Chrm(dataSize);
    std::vector<std::string> Position(dataSize);
    std::vector<std::string> Ref(dataSize);
    std::vector<std::string> Alt(dataSize);
    std::vector<double>      TaltF(dataSize);
    std::vector<int>         ReadGroupCount(dataSize);
    std::vector<bool>        IsDBSNP(dataSize);
    std::vector<std::string> RSID(dataSize);
    std::vector<bool>        UseSecondAllele(dataSize);
    std::vector<int>         TumorRefCount(dataSize);
    std::vector<int>         TumorAltCount(dataSize);
    std::vector<int>         NormalRefCount(dataSize);
    std::vector<int>         NormalAltCount(dataSize);
    std::vector<std::string> AltString(dataSize);
    std::vector<std::string> TumorFormat(dataSize);
    std::vector<std::string> NormalFormat(dataSize); 

#pragma omp parallel for
    for(int index = 0; index < dataSize; index++){
	int sorted_index = random_data_indices[index].second;
	dataEntry _dataEntry = dataParser[sorted_index];
	
	Chrm[index] = _dataEntry._Chrm;
	Position[index] = _dataEntry._Position;
	Ref[index] = _dataEntry._Ref;
	Alt[index] = _dataEntry._Alt;
	TaltF[index] = _dataEntry._TaltF;
	ReadGroupCount[index] = _dataEntry._ReadGroupCount;
	IsDBSNP[index] = "NA";
	RSID[index] = true;
	UseSecondAllele[index] = _dataEntry._UseSecondAllele;
	TumorRefCount[index] = _dataEntry._TumorRefCount; 
	TumorAltCount[index] = _dataEntry._TumorAltCount;
	NormalRefCount[index] = _dataEntry._NormalRefCount;
	NormalAltCount[index] = _dataEntry._NormalAltCount;
	AltString[index] = _dataEntry._AltString;
	TumorFormat[index] = _dataEntry._TumorFormat;
	NormalFormat[index] = _dataEntry._NormalFormat;
    } 
    
    std::vector<tabix_t *> dbsnp_vector;

    //std::cout << "num_threads:\t" << num_threads << std::endl;
    for(i = 0; i < num_threads; i++){
        //std::cout << "Start ... " << dbsnpFile << "\t" << i << std::endl;                                                                                                          
        tabix_t *t__;
        ti_iter_t iter;
        BGZF *fp;
        fp = bgzf_open(dbsnpFile, "r");
        t__ = (tabix_t *)calloc(1, sizeof(tabix_t));
        t__->fn = strdup(dbsnpFile);
        t__->fp = fp;
        if(ti_lazy_index_load(t__) < 0) {
            fprintf(stderr,"Failed to load the index of %s.\n", dbsnpFile);
            exit(EXIT_FAILURE);
        }

        dbsnp_vector.push_back(t__);
    }


#pragma omp parallel for
    for(i = 0; i < dataSize; i++){
        int _tid = omp_get_thread_num();
        tabix_t *t__ = dbsnp_vector[_tid];
        ti_iter_t iter;

        std::string rsID = "NA";
        std::string coordinate = Chrm[i] +":"+ Position[i]+"-"+Position[i];

	const char *s;
        int tid, beg, end, len;
        if(ti_parse_region(t__->idx, coordinate.c_str(), &tid, &beg, &end) == 0) {

            iter = ti_queryi(t__, tid, beg, end);

            if((s=ti_read(t__, iter, &len)) != 0) {
                std::string dbsnpInfo(s);
                std::istringstream tmpISS(dbsnpInfo);
                std::vector<std::string> tmpElements((std::istream_iterator<std::string>(tmpISS)), std::istream_iterator<std::string>());
	        rsID = tmpElements[2];
                IsDBSNP[i] = true;
                RSID[i] = rsID;
                //std::cout << i << "\t" << _tid << "\t" << rsID << std::endl;                                                                                                           
            }
            else{
                IsDBSNP[i] = false;
                RSID[i] = rsID;
            }
        }else{
            IsDBSNP[i] = false;
            RSID[i] = rsID;
        }
    }

    for(i = 0; i < num_threads; i++){
        tabix_t *t__ = dbsnp_vector[i];
        ti_close(t__);
    }

    
    // check if the BAM has multiple read group
    //
    bool isMultiReadGroup = false;
    for(i = 0; i < ReadGroupCount.size(); i++) {
        if(ReadGroupCount[i] > 1) {
            isMultiReadGroup = true;
            break;
        }
    }

    // process "MuSE call" header information
    //
    std::vector<std::string> MuSEVersion;
    std::vector<std::string> MuSECallCMD;
    std::vector<std::string> TumorSample;
    std::vector<std::string> NormalSample;
    std::vector<std::string> ContigInfo;
    std::vector<std::string> ReferenceGenome;
    for(i = 0; i < Header.size(); i++) {
        if(Header[i].find("Build Date")!=std::string::npos && Header[i].find("Build Time")!=std::string::npos) {
            if(MuSEVersion.empty()) {
                MuSEVersion.push_back(Header[i]);
            }
            else {
                if(std::find(MuSEVersion.begin(), MuSEVersion.end(), Header[i])==MuSEVersion.end()) {
                    MuSEVersion.push_back(Header[i]);
                }
            }
        }
        else if(Header[i].find("##MuSE_call")!=std::string::npos) {
            if(MuSECallCMD.empty()) {
                MuSECallCMD.push_back(Header[i]);
            }
            else {
                if(std::find(MuSECallCMD.begin(), MuSECallCMD.end(), Header[i])==MuSECallCMD.end()) {
                    MuSECallCMD.push_back(Header[i]);
                }
            }
        }
        else if(Header[i].find("##TUMOR")!=std::string::npos) {
            if(TumorSample.empty()) {
                TumorSample.push_back(Header[i]);
            }
            else {
                if(std::find(TumorSample.begin(), TumorSample.end(), Header[i])==TumorSample.end()) {
                    TumorSample.push_back(Header[i]);
                }
            }
        }
        else if(Header[i].find("##NORMAL")!=std::string::npos) {
            if(NormalSample.empty()) {
                NormalSample.push_back(Header[i]);
            }
            else {
                if(std::find(NormalSample.begin(), NormalSample.end(), Header[i])==NormalSample.end()) {
                    NormalSample.push_back(Header[i]);
                }
            }
        }
        else if(Header[i].find("##contig=<ID=")!=std::string::npos) {
            if(ContigInfo.empty()) {
                ContigInfo.push_back(Header[i]);
            }
            else {
                if(std::find(ContigInfo.begin(), ContigInfo.end(), Header[i])==ContigInfo.end()) {
                    ContigInfo.push_back(Header[i]);
                }
            }
        }
        else if(Header[i].find("##reference=file")!=std::string::npos) {
            if(ReferenceGenome.empty()) {
                ReferenceGenome.push_back(Header[i]);
            }
            else {
                if(std::find(ReferenceGenome.begin(), ReferenceGenome.end(), Header[i])==ReferenceGenome.end()) {
                    ReferenceGenome.push_back(Header[i]);
                }
            }
        }
    }

    // check extracted header info
    //
    if(MuSEVersion.size()>1) {
        std::cerr << "Multiple MuSE versions found. Please check input file.\n";
        exit(EXIT_FAILURE);
    }
    if(TumorSample.size()!=1) {
        std::cerr << "Either tumor sample ID is missing or multiple tumor sample IDs are found. Please check input file.\n";
        exit(EXIT_FAILURE);
    }
    if(NormalSample.size()!=1) {
        std::cerr << "Either normal sample ID is missing or multiple normal sample IDs are found. Please check input file.\n";
        exit(EXIT_FAILURE);
    }
    if(ReferenceGenome.size()>1) {
        std::cerr << "Multiple reference genomes used. Please check input file.\n";
        exit(EXIT_FAILURE);
    }

    // extract contig names, used for sorting "MuSE call" output
    //

    std::vector<std::string> ContigName;
    for(i = 0; i < ContigInfo.size(); i++) {
        std::istringstream ss(ContigInfo[i]);
        std::string token;
        std::getline(ss, token, ',');
        ContigName.push_back(token.substr(13));
    }
    
    // sort extracted "MuSE call" output by coordinate
    //
    std::vector<std::string> FinalChrm;
    std::vector<int>         FinalPosition;
    std::vector<std::string> FinalRef;
    std::vector<std::string> FinalAlt;
    std::vector<double>      FinalTaltF;
    std::vector<int>         FinalReadGroupCount;
    std::vector<bool>        FinalIsDBSNP;
    std::vector<std::string> FinalRSID;
    std::vector<bool>        FinalUseSecondAllele;
    std::vector<int>         FinalTumorRefCount;
    std::vector<int>         FinalTumorAltCount;
    std::vector<int>         FinalNormalRefCount;
    std::vector<int>         FinalNormalAltCount;
    std::vector<std::string> FinalAltString;
    std::vector<std::string> FinalTumorFormat;
    std::vector<std::string> FinalNormalFormat;

    FinalChrm.reserve(dataSize);
    FinalPosition.reserve(dataSize);
    FinalRef.reserve(dataSize);
    FinalAlt.reserve(dataSize);
    FinalTaltF.reserve(dataSize);
    FinalReadGroupCount.reserve(dataSize);
    FinalIsDBSNP.reserve(dataSize);
    FinalRSID.reserve(dataSize);
    FinalUseSecondAllele.reserve(dataSize);
    FinalTumorRefCount.reserve(dataSize);
    FinalTumorAltCount.reserve(dataSize);
    FinalNormalRefCount.reserve(dataSize);
    FinalNormalAltCount.reserve(dataSize);
    FinalAltString.reserve(dataSize);
    FinalTumorFormat.reserve(dataSize);
    FinalNormalFormat.reserve(dataSize);

    typedef std::map<int, callInfo>::value_type callType;
    
    for(i = 0; i < ContigName.size(); i++) {
        std::map<int, callInfo, std::less<int>> ContigCall;
        //std::vector<int> TmpPosition;
        //TmpPosition.reserve(100000);

#pragma omp parallel for 
        for(j = 0; j < Chrm.size(); j++) {
            if(Chrm[j]==ContigName[i]) {
                int position = atoi(Position[j].c_str());

                callInfo oneCall;
                oneCall.ref             = Ref[j];
                oneCall.alt             = Alt[j];
                oneCall.tAltF           = TaltF[j];
                oneCall.readGroupCount  = ReadGroupCount[j];
                oneCall.isDBSNP         = IsDBSNP[j];
                oneCall.rsID            = RSID[j];
                oneCall.useSecondAllele = UseSecondAllele[j];
                oneCall.tumorRefCount   = TumorRefCount[j];
                oneCall.tumorAltCount   = TumorAltCount[j];
                oneCall.normalRefCount  = NormalRefCount[j];
                oneCall.normalAltCount  = NormalAltCount[j];
                oneCall.altString       = AltString[j];
                oneCall.tumorFormat     = TumorFormat[j];
                oneCall.normalFormat    = NormalFormat[j];

#pragma omp critical
		{
		    //TmpPosition.push_back(position);
		    ContigCall.insert(callType(position, oneCall));
		}
            }
        }
	
        //std::sort(TmpPosition.begin(), TmpPosition.end());

	for(auto &_callType: ContigCall){

	    callInfo oneCall = _callType.second;
	    //for(j = 0; j < TmpPosition.size(); j++) {
            FinalChrm.push_back(ContigName[i]);
            FinalPosition.push_back(_callType.first);
            FinalRef.push_back(oneCall.ref);
            FinalAlt.push_back(oneCall.alt);
            FinalTaltF.push_back(oneCall.tAltF);
            FinalReadGroupCount.push_back(oneCall.readGroupCount);
            FinalIsDBSNP.push_back(oneCall.isDBSNP);
            FinalRSID.push_back(oneCall.rsID);
            FinalUseSecondAllele.push_back(oneCall.useSecondAllele);
            FinalTumorRefCount.push_back(oneCall.tumorRefCount);
            FinalTumorAltCount.push_back(oneCall.tumorAltCount);
            FinalNormalRefCount.push_back(oneCall.normalRefCount);
            FinalNormalAltCount.push_back(oneCall.normalAltCount);
            FinalAltString.push_back(oneCall.altString);
            FinalTumorFormat.push_back(oneCall.tumorFormat);
            FinalNormalFormat.push_back(oneCall.normalFormat);
	    //}
	}
    }

    // identify the dynamic cutoff, either Beta or GMM based on the input option
    //
    

    double tier1Cutoff = 10.0;
    double tier2Cutoff = 10.0;
    double tier3Cutoff = 10.0;
    double tier4Cutoff = 10.0;

    bool isWGSFailed = false;

    if(isWGS) {
        // prepare the data
        //
        std::vector<double> lnTEVAF;

        lnTEVAF.reserve(dataSize);
        if(isMultiReadGroup) {
            for(i = 0; i < dataSize; i++) {
                if(FinalReadGroupCount[i] > 1) {
                    lnTEVAF.push_back(log(FinalTaltF[i]));
                }
            }
        }
        else {
            for(i = 0; i < dataSize; i++) {
                lnTEVAF.push_back(log(FinalTaltF[i]));
            }
        }
        std::sort(lnTEVAF.begin(), lnTEVAF.end());

        // EM based on mixtools v1.0.2 normalmixEM, always assume two components
        //
        int    emReplicate  = 50;
        int    maxIteration = 5000;
        int    maxRestarts  = 20;
        int    dataSize     = (int)lnTEVAF.size();
        double epsilon      = 1e-08;

        // partition
        //
        int maxIndexLeftPartition  = (int)ceil((double)dataSize/2.0) - 1;
        int minIndexRightPartition = (int)floor((double)dataSize/2.0) - 1;
        // mean
        //
        double leftPartitionMean = 0.0;
        for(i = 0; i <= maxIndexLeftPartition; i++) {
            leftPartitionMean += lnTEVAF[i];
        }
        leftPartitionMean = leftPartitionMean/(double)(maxIndexLeftPartition+1);
        double rightPartitionMean = 0.0;
        for(i = minIndexRightPartition; i < dataSize; i++) {
            rightPartitionMean += lnTEVAF[i];
        }
        rightPartitionMean = rightPartitionMean/(double)(dataSize-minIndexRightPartition);
        double lnTEVAFMean = 0.0;
        for(i = 0; i < dataSize; i++) {
            lnTEVAFMean += lnTEVAF[i];
        }
        lnTEVAFMean = lnTEVAFMean/(double)dataSize;
        // standard deviation
        //
        double leftPartitionSD = 0.0;

#pragma omp parallel for reduction(+:leftPartitionSD) 
        for(i = 0; i <= maxIndexLeftPartition; i++) {
            leftPartitionSD += (lnTEVAF[i] - leftPartitionMean)*(lnTEVAF[i] - leftPartitionMean);
        }
        leftPartitionSD = sqrt(leftPartitionSD/(double)(maxIndexLeftPartition));
        double originalLeftPartitionSD = leftPartitionSD;
        double rightPartitionSD = 0.0;
	
#pragma omp parallel for reduction(+:rightPartitionSD) 
        for(i = minIndexRightPartition; i < dataSize; i++) {
            rightPartitionSD += (lnTEVAF[i] - rightPartitionMean)*(lnTEVAF[i] - rightPartitionMean);
        }
        rightPartitionSD = sqrt(rightPartitionSD/(double)(dataSize-minIndexRightPartition-1));
        double originalRightPartitionSD = rightPartitionSD;
        double lnTEVAFSD = 0.0;

#pragma omp parallel for reduction(+:lnTEVAFSD) 
        for(i = 0; i < dataSize; i++) {
            lnTEVAFSD += (lnTEVAF[i] - lnTEVAFMean)*(lnTEVAF[i] - lnTEVAFMean);
        }
	
        lnTEVAFSD = sqrt(lnTEVAFSD/(double)(dataSize-1));

	        // try EM 50 times
        //

	//std::cout << "Original:" << "\t" << leftPartitionSD << "\t" << rightPartitionSD << std::endl;	    
        EMEstimate.clear();

#pragma omp parallel for private(i, j)
        for(ithReplicate = 0; ithReplicate < emReplicate; ithReplicate++) {
            // initialize mu
            //
            double mu[2] = {-6.0, -2.0};
            // propose sigma
            //
            double sigma[2];
            if(originalLeftPartitionSD==0.0) {
                do {leftPartitionSD = unif_rand()*lnTEVAFSD;} while(leftPartitionSD == 0.0);
            }
            if(originalRightPartitionSD==0.0) {
                do {rightPartitionSD = unif_rand()*lnTEVAFSD;} while(rightPartitionSD == 0.0);
            }
            sigma[0] = leftPartitionSD/exp_rand();
            sigma[1] = rightPartitionSD/exp_rand();
            //propose lambda
            //

            double lambda[2];
	    double tmpLambda1 = unif_rand();
	    double tmpLambda2 = unif_rand();
	    lambda[0] = tmpLambda1/(tmpLambda1 + tmpLambda2);
	    lambda[1] = tmpLambda2/(tmpLambda1 + tmpLambda2);
	    // monitor if done and count the number of restart
	    // 
	    
            bool   notDone   = true;
            int    nRestart  = 0;
            int    iteration = 0;
            double loglik    = 0.0;
            while(notDone) {
                // initialize
                //
                notDone = false;
		
                double diff = epsilon + 1;
                iteration = 0;
                std::vector<double> postProbs(dataSize*2, 0.0);
                std::vector<double> res(dataSize*2, 0.0);
                std::vector<double> work(3*2, 0.0);
                loglik = 0.0;
                // Initialization E-step
                //
                normpost(dataSize, &lnTEVAF[0], mu, sigma, lambda, &res[0], &work[0], &postProbs[0], &loglik);
                double obsloglik = loglik;
                while(diff > epsilon && iteration < maxIteration) {
                    // M-step
                    //
                    for(i = 0; i < 2; i++) {
                        lambda[i] = 0.0;
                        for(j = 0; j < dataSize; j++) {
                            lambda[i] += postProbs[i*dataSize+j];
                        }
                        lambda[i] = lambda[i]/(double)dataSize;
                    }
                    for(i = 0; i < 2; i++) {
                        mu[i] = 0.0;
                        for(j = 0; j < dataSize; j++) {
                            mu[i] += postProbs[i*dataSize+j] * lnTEVAF[j];
                        }
                        mu[i] = mu[i]/((double)dataSize*lambda[i]);
                    }
                    for(i = 0; i < 2; i++) {
                        sigma[i] = 0.0;
                        for(j = 0; j < dataSize; j++) {
                            sigma[i] += postProbs[i*dataSize+j] * res[i*dataSize+j];
                        }
                        sigma[i] = sqrt(sigma[i]/((double)dataSize*lambda[i]));
                    }
                    if(sigma[0]<1e-08 || sigma[1]<1e-08) {
                        notDone = true;
                        std::cout << "One of the variances is going to zero; trying new starting values.\n";
                        nRestart += 1;
                        // propose new lambda, mu and sigma
                        //
                        tmpLambda1 = unif_rand();
                        tmpLambda2 = unif_rand();
                        lambda[0]  = tmpLambda1/(tmpLambda1 + tmpLambda2);
                        lambda[1]  = tmpLambda2/(tmpLambda1 + tmpLambda2);
                        if(originalLeftPartitionSD==0.0) {
                            do {leftPartitionSD = unif_rand()*lnTEVAFSD;} while(leftPartitionSD == 0.0);
                        }
                        if(originalRightPartitionSD==0.0) {
                            do {rightPartitionSD = unif_rand()*lnTEVAFSD;} while(rightPartitionSD == 0.0);
                        }
                        sigma[0] = leftPartitionSD/exp_rand();
                        sigma[1] = rightPartitionSD/exp_rand();
                        mu[0] = leftPartitionMean  + sigma[0]*norm_rand();
                        mu[1] = rightPartitionMean + sigma[1]*norm_rand();
                        if(nRestart>maxRestarts) {
                            std::cout << "Too many tries!\n";
                            exit(EXIT_FAILURE);
                        }
                        break;
                    }
                    // E-step
                    //
		    //std::cout << "iteration: " << iteration << "\t" << "loglik: " << loglik << std::endl;
                    normpost(dataSize, &lnTEVAF[0], mu, sigma, lambda, &res[0], &work[0], &postProbs[0], &loglik);
                    double newobsloglik = loglik;
                    diff = newobsloglik - obsloglik;
                    obsloglik = newobsloglik;
                    iteration += 1;
                }
            }
            if(iteration==maxIteration) {
                std::cout << "WARNING! NOT CONVERGENT!\n";
            }

	    // mu[0] and mu[1] cannot be equal
            //
            //if(mu[0]==mu[1])
	    //   continue;
            // make sure mu[0], sigma[0] and lambda[0] associate with left normal
            //
            if(mu[0]>mu[1]) {
                double tmp = mu[0];
                mu[0]     = mu[1];
                mu[1]     = tmp;
                tmp       = sigma[0];
                sigma[0]  = sigma[1];
                sigma[1]  = tmp;
                tmp       = lambda[0];
                lambda[0] = lambda[1];
                lambda[1] = tmp;
            }
            // check the EM estimates
            // put some constraints here, if meet, then record it
            
	    
            if(mu[0]!=mu[1] &&  mu[1]>-4.605170185988091) {
                emEst tmpEstimate;
                tmpEstimate.loglik  = loglik;
                tmpEstimate.mu1     = mu[0];
                tmpEstimate.mu2     = mu[1];
                tmpEstimate.sigma1  = sigma[0];
                tmpEstimate.sigma2  = sigma[1];
                tmpEstimate.lambda1 = lambda[0];
                tmpEstimate.lambda2 = lambda[1];

#pragma omp critical
		{
		    EMEstimate.push_back(tmpEstimate);
		}
            }
        }

	std::cout <<"EMEstimate size: " << EMEstimate.size() << std::endl;

        // calculate the mean mu from 50 replicates
        // then, fix mu to re-estimate the GMM
        //
        if(EMEstimate.size()==0) {
            std::cout << "Not enough data points for model fitting. Automatically switch to option -E.\n";
            isWGSFailed = true;
        }
        else {
            // constrained mu
            //
            finalMu[0] = 0.0;
            finalMu[1] = 0.0;
            for(i = 0; i < EMEstimate.size(); i++) {
                finalMu[0] += EMEstimate[i].mu1;
                finalMu[1] += EMEstimate[i].mu2;
            }
            finalMu[0] = finalMu[0]/(double)EMEstimate.size();
            finalMu[1] = finalMu[1]/(double)EMEstimate.size();
            // propose sigma
            //
            if(originalLeftPartitionSD==0.0) {
                do {leftPartitionSD = unif_rand()*lnTEVAFSD;} while(leftPartitionSD == 0.0);
            }
            if(originalRightPartitionSD==0.0) {
                do {rightPartitionSD = unif_rand()*lnTEVAFSD;} while(rightPartitionSD == 0.0);
            }
            finalSigma[0] = leftPartitionSD/exp_rand();
            finalSigma[1] = rightPartitionSD/exp_rand();
            //propose lambda
            //
            double tmpLambda1 = unif_rand();
            double tmpLambda2 = unif_rand();
            finalLambda[0] = tmpLambda1/(tmpLambda1 + tmpLambda2);
            finalLambda[1] = tmpLambda2/(tmpLambda1 + tmpLambda2);
            // monitor if done and count the number of restart
            //
            bool   notDone   = true;
            int    nRestart  = 0;
            int    iteration = 0;
            double loglik    = 0.0;
            while(notDone) {
                // initialize
                //
                notDone = false;
                double diff = epsilon + 1;
                iteration = 0;
                std::vector<double> postProbs(dataSize*2, 0.0);
                std::vector<double> res(dataSize*2, 0.0);
                std::vector<double> work(3*2, 0.0);
                loglik = 0.0;
                // Initialization E-step
                //
                normpost(dataSize, &lnTEVAF[0], finalMu, finalSigma, finalLambda, &res[0], &work[0], &postProbs[0], &loglik);
                double obsloglik = loglik;
                while(diff > epsilon && iteration < maxIteration) {
                    // M-step
                    //
                    for(i = 0; i < 2; i++) {
                        finalLambda[i] = 0.0;
                        for(j = 0; j < dataSize; j++) {
                            finalLambda[i] += postProbs[i*dataSize+j];
                        }
                        finalLambda[i] = finalLambda[i]/(double)dataSize;
                    }
                    for(i = 0; i < 2; i++) {
                        finalSigma[i] = 0.0;
                        for(j = 0; j < dataSize; j++) {
                            finalSigma[i] += postProbs[i*dataSize+j] * res[i*dataSize+j];
                        }
                        finalSigma[i] = sqrt(finalSigma[i]/((double)dataSize*finalLambda[i]));
                    }
                    if(finalSigma[0]<1e-08 || finalSigma[1]<1e-08) {
                        notDone = true;
                        std::cout << "One of the variances is going to zero; trying new starting values.\n";
                        nRestart += 1;
                        // propose new lambda and sigma
                        //
                        tmpLambda1 = unif_rand();
                        tmpLambda2 = unif_rand();
                        finalLambda[0] = tmpLambda1/(tmpLambda1 + tmpLambda2);
                        finalLambda[1] = tmpLambda2/(tmpLambda1 + tmpLambda2);
                        if(originalLeftPartitionSD==0.0) {
                            do {leftPartitionSD = unif_rand()*lnTEVAFSD;} while(leftPartitionSD == 0.0);
                        }
                        if(originalRightPartitionSD==0.0) {
                            do {rightPartitionSD = unif_rand()*lnTEVAFSD;} while(rightPartitionSD == 0.0);
                        }
                        finalSigma[0] = leftPartitionSD/exp_rand();
                        finalSigma[1] = rightPartitionSD/exp_rand();
                        if(nRestart>maxRestarts) {
                            std::cout << "Too many tries!\n";
                            exit(EXIT_FAILURE);
                        }
                        break;
                    }
                    // E-step
                    //
                    normpost(dataSize, &lnTEVAF[0], finalMu, finalSigma, finalLambda, &res[0], &work[0], &postProbs[0], &loglik);
                    double newobsloglik = loglik;
                    diff = newobsloglik - obsloglik;
                    obsloglik = newobsloglik;
                    iteration += 1;
                }
            }
            if(iteration==maxIteration) {
                std::cout << "WARNING! NOT CONVERGENT!\n";
            }
            // assume EM finished successfully, use mu, sigma, lambda to find the cutoff that minimizes the misclassification prob
            //
            int    n_round     = 0;
            int    n_round_max = 500;
            double prob_new    = UNLIKELY;
            double prob_old    = UNLIKELY;
            gmmCutoff = 0.5*(finalMu[0] + finalMu[1]);
            Negative_Misclassification_Prob();
            while(n_round < n_round_max) {
                Generic_Brent_Lk(&gmmCutoff,
                                 finalMu[0],
                        finalMu[1],
                        min_diff_misclassification_prob,
                        BRENT_IT_MAX,
                        0,
                        Negative_Misclassification_Prob);
                prob_new = negative_misclassification_prob;
                if(fabs(prob_new - prob_old) < min_diff_misclassification_prob)
                    break;
                else
                    prob_old = prob_new;
                n_round++;
            }

            // comapre exp(gmmCutoff) with passVAF 0.02, the larger one will be the new passVAF, and the smaller one will be the tier1 cutoff
            //
            if(exp(gmmCutoff) >= passVAF) {
                // the top01 here is the tier1 cutoff
                //
                tier1Cutoff = passVAF;
                passVAF = exp(gmmCutoff);
            }
            else {
                tier1Cutoff = exp(gmmCutoff);
            }

            // the tier2, tier3 and tier4 cutoffs will be the top 0.1, 0.5 and 1 percentile of the truncated normal (truncation ends at the tier1 cutoff)
            // area under the left Gaussian curve up to the tier1 cutoff
            //
            double leftArea  = pnorm(log(tier1Cutoff), finalMu[0], finalSigma[0], 1, 0);
            double rightArea = 1.0 - leftArea;
            tier2Cutoff = exp(qnorm(rightArea+0.001*leftArea, finalMu[0], finalSigma[0], 0, 0));
            tier3Cutoff = exp(qnorm(rightArea+0.005*leftArea, finalMu[0], finalSigma[0], 0, 0));
            tier4Cutoff = exp(qnorm(rightArea+0.010*leftArea, finalMu[0], finalSigma[0], 0, 0));
        }
    } 
	if(isWES || (isWGS && isWGSFailed)) {
		// prepare the input data
		//
		altf.clear();
		altf.reserve(dataSize);
		for(i = 0; i < dataSize; i++) {
			if(FinalTaltF[i]<=passVAF && FinalTaltF[i]>=minTumorVAF) {
				double tmp = (FinalTaltF[i] - minTumorVAF)/(passVAF - minTumorVAF);
				if(tmp > 0.0 && tmp < 1.0) {
					altf.push_back(tmp);
				}
			}
		}
		
		// mle of beta shape parameters
		//
		int    n_round     = 0;
		int    n_round_max = 500;
		double lk_new      = UNLIKELY;
		double lk_old      = UNLIKELY;
		Beta_Lk();
		while(n_round < n_round_max) {
		  //int i;
			For(i,2) {
				Generic_Brent_Lk(&(betaShape[i]),
								 0.,
								 100.,
								 min_diff_beta_lk,
								 BRENT_IT_MAX,
								 0,
								 Beta_Lk);
			}
			lk_new = beta_lnlike;
			if(fabs(lk_new - lk_old) < min_diff_beta_lk)
				break;
			else
				lk_old = lk_new;
			n_round++;
		}
		
		// tier-based cutoff
		//
		tier1Cutoff = qbeta(0.999, betaShape[0], betaShape[1], 1, 0)*(passVAF-minTumorVAF)+minTumorVAF;
		tier2Cutoff = qbeta(0.995, betaShape[0], betaShape[1], 1, 0)*(passVAF-minTumorVAF)+minTumorVAF;
		tier3Cutoff = qbeta(0.990, betaShape[0], betaShape[1], 1, 0)*(passVAF-minTumorVAF)+minTumorVAF;
		tier4Cutoff = qbeta(0.980, betaShape[0], betaShape[1], 1, 0)*(passVAF-minTumorVAF)+minTumorVAF;
	}
	
#if DEBUG
	if(isWGS && (!isWGSFailed)) {
		std::cout << "successful tries:\n";
		for(i = 0; i < EMEstimate.size(); i++) {
			std::cout << EMEstimate[i].loglik << '\t' << EMEstimate[i].mu1 << '\t' << EMEstimate[i].mu2 << '\t' << EMEstimate[i].sigma1 << '\t' << EMEstimate[i].sigma2 << '\t' << EMEstimate[i].lambda1 << '\t' << EMEstimate[i].lambda2 << '\n';
		}
		std::cout << "mu : " << finalMu[0] << '\t' << finalMu[1] << '\n';
		std::cout << "sigma : " << finalSigma[0] << '\t' << finalSigma[1] << '\n';
		std::cout << "lambda : " << finalLambda[0] << '\t' << finalLambda[1] << '\n';
		std::cout << "Dynamic Cutoff : " << exp(gmmCutoff) << '\t' << tier2Cutoff << '\t' << tier3Cutoff << '\t' << tier4Cutoff << '\n';
	}
	if(isWES || (isWGS && isWGSFailed)) {
		std::cout << "Beta Shape : " << betaShape[0] << '\t' << betaShape[1] << '\n';
		std::cout << "Dynamic Cutoff : " << tier1Cutoff << '\t' << tier2Cutoff << '\t' << tier3Cutoff << '\t' << tier4Cutoff << '\n';
	}
#endif

	// generate final call set in VCF format
	//
	// write out to VCF
	//
	std::ofstream outf(outFile);
	outf << "##fileformat=VCFv4.1\n";
	outf << "##FILTER=<ID=PASS,Description=\"Accept as a confident somatic mutation\">\n";
	outf << "##FILTER=<ID=Tier1,Description=\"Confident level 1\">\n";
	outf << "##FILTER=<ID=Tier2,Description=\"Confident level 2\">\n";
	outf << "##FILTER=<ID=Tier3,Description=\"Confident level 3\">\n";
	outf << "##FILTER=<ID=Tier4,Description=\"Confident level 4\">\n";
	outf << "##FILTER=<ID=Tier5,Description=\"Confident level 5\">\n";
	outf << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	outf << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth at this position in the sample\">\n";
	outf << "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Depth of reads supporting alleles 0/1/2/3...\">\n";
	outf << "##FORMAT=<ID=BQ,Number=.,Type=Integer,Description=\"Average base quality for reads supporting alleles\">\n";
	outf << "##FORMAT=<ID=SS,Number=1,Type=Integer,Description=\"Variant status relative to non-adjacent Normal,0=wildtype,1=germline,2=somatic,3=LOH,4=post-transcriptional modification,5=unknown\">\n";
	outf << "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Indicates if record is a somatic mutation\">\n";
	outf << TumorSample[0] << '\n';
	outf << NormalSample[0] << '\n';
	for(i = 0; i < MuSEVersion.size(); i++) {
		outf << MuSEVersion[i] << '\n';
	}
	for(i = 0; i < MuSECallCMD.size(); i++) {
		std::string base = MuSECallCMD[i];
		std::string substr = "##MuSE_call_" + IntToString(i+1);
		base.replace(0, 11, substr);
		outf << base << '\n';
	}
	// sump command line
	//
	outf << "##MuSE_sump=\"";
	for(i = 0; i < argc-1; i++) {
		outf << argv[i] << " ";
	}
	outf << argv[argc-1] << "\"\n";
	for(i = 0; i < ContigInfo.size(); i++) {
		outf << ContigInfo[i] << '\n';
	}
	for(i = 0; i < ReferenceGenome.size(); i++) {
		outf << ReferenceGenome[i] << '\n';
	}
	outf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL\n";
	std::string whichTier = "";
	for(i = 0; i < FinalTaltF.size(); i++) {
		if(FinalIsDBSNP[i]) {
			if(FinalTaltF[i] >= 2.0*passVAF) {
				whichTier = "PASS";
				if(isMultiReadGroup) {
					if(FinalReadGroupCount[i]>1)
						outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << FinalRSID[i] << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
				else {
					outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << FinalRSID[i] << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
			}
			else if(FinalTaltF[i] >= 2.0*tier1Cutoff && FinalTaltF[i] >= 2.0*baseVAF) {
				if(tier1Cutoff > baseVAF) {
					whichTier = "Tier1";
				}
				else {
					whichTier = "Tier5";
				}
				if(isMultiReadGroup) {
					if(FinalReadGroupCount[i]>1)
						outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << FinalRSID[i] << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
				else {
					outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << FinalRSID[i] << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
			}
			else if(FinalTaltF[i] >= 2.0*tier2Cutoff && FinalTaltF[i] >= 2.0*baseVAF) {
				if(tier2Cutoff > baseVAF) {
					whichTier = "Tier2";
				}
				else {
					whichTier = "Tier5";
				}
				if(isMultiReadGroup) {
					if(FinalReadGroupCount[i]>1)
						outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << FinalRSID[i] << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
				else {
					outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << FinalRSID[i] << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
			}
			else if(FinalTaltF[i] >= 2.0*tier3Cutoff && FinalTaltF[i] >= 2.0*baseVAF) {
				if(tier3Cutoff > baseVAF) {
					whichTier = "Tier3";
				}
				else {
					whichTier = "Tier5";
				}
				if(isMultiReadGroup) {
					if(FinalReadGroupCount[i]>1)
						outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << FinalRSID[i] << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
				else {
					outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << FinalRSID[i] << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
			}
			else if(FinalTaltF[i] >= 2.0*tier4Cutoff && FinalTaltF[i] >= 2.0*baseVAF) {
				if(tier4Cutoff > baseVAF) {
					whichTier = "Tier4";
				}
				else {
					whichTier = "Tier5";
				}
				if(isMultiReadGroup) {
					if(FinalReadGroupCount[i]>1)
						outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << FinalRSID[i] << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
				else {
					outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << FinalRSID[i] << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
			}
			else if(FinalTaltF[i] >= 2.0*baseVAF) {
				whichTier = "Tier5";
				if(isMultiReadGroup) {
					if(FinalReadGroupCount[i]>1)
						outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << FinalRSID[i] << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
				else {
					outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << FinalRSID[i] << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
			}
		}
		else {
			if(FinalTaltF[i] >= passVAF) {
				whichTier = "PASS";
				if(isMultiReadGroup) {
					if(FinalReadGroupCount[i]>1)
						outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << "." << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
				else {
					outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << "." << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
			}
			else if(FinalTaltF[i] >= tier1Cutoff && FinalTaltF[i] >= baseVAF) {
				if(tier1Cutoff > baseVAF) {
					whichTier = "Tier1";
				}
				else {
					whichTier = "Tier5";
				}
				if(isMultiReadGroup) {
					if(FinalReadGroupCount[i]>1)
						outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << "." << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
				else {
					outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << "." << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
			}
			else if(FinalTaltF[i] >= tier2Cutoff && FinalTaltF[i] >= baseVAF) {
				if(tier2Cutoff > baseVAF) {
					whichTier = "Tier2";
				}
				else {
					whichTier = "Tier5";
				}
				if(isMultiReadGroup) {
					if(FinalReadGroupCount[i]>1)
						outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << "." << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
				else {
					outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << "." << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
			}
			else if(FinalTaltF[i] >= tier3Cutoff && FinalTaltF[i] >= baseVAF) {
				if(tier3Cutoff > baseVAF) {
					whichTier = "Tier3";
				}
				else {
					whichTier = "Tier5";
				}
				if(isMultiReadGroup) {
					if(FinalReadGroupCount[i]>1)
						outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << "." << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
				else {
					outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << "." << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
			}
			else if(FinalTaltF[i] >= tier4Cutoff && FinalTaltF[i] >= baseVAF) {
				if(tier4Cutoff > baseVAF) {
					whichTier = "Tier4";
				}
				else {
					whichTier = "Tier5";
				}
				if(isMultiReadGroup) {
					if(FinalReadGroupCount[i]>1)
						outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << "." << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
				else {
					outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << "." << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
			}
			else if(FinalTaltF[i] >= baseVAF) {
				whichTier = "Tier5";
				if(isMultiReadGroup) {
					if(FinalReadGroupCount[i]>1)
						outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << "." << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
				else {
					outf << FinalChrm[i] << '\t' << FinalPosition[i] << '\t' << "." << '\t' << FinalRef[i] << '\t' << FinalAltString[i] << '\t' << "." << '\t' << whichTier << "\tSOMATIC\tGT:DP:AD:BQ:SS\t" << FinalTumorFormat[i] << '\t' << FinalNormalFormat[i] << '\n';
				}
			}
		}
	}
	
	outf.close();
}
