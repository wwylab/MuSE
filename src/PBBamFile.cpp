

#include <assert.h>
#include <string>
#include <vector>
#include <iostream>
#include "PBBamFile.h"
#include "htslib/khash.h"
#include "htslib/bgzf.h"


//READING FILE
/* The following macro calls a zlib routine and checks the return
   value. If the return value ("status") is not OK, it prints an error
   message and exits the program. Zlib's error statuses are all less
   than zero. */

#define CALL_ZLIB(x) {                                                  \
        int status;                                                     \
        status = x;                                                     \
        if (status < 0) {                                               \
            fprintf (stderr,                                            \
                     "%s:%d: %s returned a bad status of %d.\n",        \
                     __FILE__, __LINE__, #x, status);                   \
            exit (EXIT_FAILURE);                                        \
        }                                                               \
    }
#define CHUNKSIZE 0xff00
#define READCHUNKSIZE 0x10000
#define MAXCOMPRESSEDDATASIZE 100000000

using namespace std;

 void swap_data(const bam1_core_t *c, int l_data, uint8_t *data, int is_host)
{
    uint32_t *cigar = (uint32_t*)(data + c->l_qname);
    uint32_t i;
    for (i = 0; i < c->n_cigar; ++i) ed_swap_4p((void*)&cigar[i]);
}

uint32_t fromByteArray(uint8_t * bytes) {
    uint32_t result = (bytes[3]<<24) | (bytes[2]<<16) | (bytes[1]<<8) | bytes[0];
    return result;
}
 
 inline uint32_t pb_ed_swap_4(uint32_t v)
{
    v = ((v & 0x0000FFFFU) << 16) | (v >> 16);
    return ((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8);
}
 inline void *pb_ed_swap_4p(void *x)
{
    *(uint32_t*)x = pb_ed_swap_4(*(uint32_t*)x);
    return x;
}

//#define pair64_lt(a,b) ((a).u < (b).u)

//KSORT_INIT(_off, hts_pair64_t, pair64_lt)

typedef struct {
    int32_t m, n;
    uint64_t loff;
    hts_pair64_t *list;
} bins_t;

KHASH_MAP_INIT_INT(bin, bins_t)
typedef khash_t(bin) bidx_t;

typedef struct {
    int32_t n, m;
    uint64_t *offset;
} lidx_t;

struct __hts_idx_t {
    int fmt, min_shift, n_lvls, n_bins;
    uint32_t l_meta;
    int32_t n, m;
    uint64_t n_no_coor;
    bidx_t **bidx;
    lidx_t *lidx;
    uint8_t *meta; // MUST have a terminating NUL on the end
    struct {
        uint32_t last_bin, save_bin;
        int last_coor, last_tid, save_tid, finished;
        uint64_t last_off, save_off;
        uint64_t off_beg, off_end;
        uint64_t n_mapped, n_unmapped;
    } z; // keep internal states
};

int bam_read1_pb(bam1_t *b, uint8_t * data, uint32_t size, bool is_be)
{
    bam1_core_t *c = &b->core;
    int cur = 0;
    int32_t block_len, ret, i;
    uint32_t x[8];
    if (cur + 4 > size)
        return -2;
    block_len = fromByteArray(data + cur);
    cur += 4;

    if (is_be) {
        pb_ed_swap_4p(&block_len);
    }

    if (cur + block_len > size)
        return -3;
    for (int i = 0; i < 8; i++)  {  
        x[i] = fromByteArray(data + cur);
        cur += 4;
    }
    
    if (is_be) {     
        for (i = 0; i < 8; ++i) pb_ed_swap_4p(x + i);
    }
    c->tid = x[0]; c->pos = x[1];
    c->bin = x[2]>>16; c->qual = x[2]>>8&0xff; c->l_qname = x[2]&0xff;
    c->l_extranul = (c->l_qname%4 != 0)? (4 - c->l_qname%4) : 0;
    if ((uint32_t) c->l_qname + c->l_extranul > 255) // l_qname would overflow
        return -4;
    c->flag = x[3]>>16; c->n_cigar = x[3]&0xffff;
    c->l_qseq = x[4];
    c->mtid = x[5]; c->mpos = x[6]; c->isize = x[7];
    b->l_data = block_len - 32 + c->l_extranul;
    if (b->l_data < 0 || c->l_qseq < 0 || c->l_qname < 1) return -4;
    if (((uint64_t) c->n_cigar << 2) + c->l_qname + c->l_extranul
        + (((uint64_t) c->l_qseq + 1) >> 1) + c->l_qseq > (uint64_t) b->l_data)
        return -4;    
    if (b->m_data < b->l_data) {
        uint8_t *new_data;
        uint32_t new_m = b->l_data;
        kroundup32(new_m);
        new_data = (uint8_t*)realloc(b->data, new_m);
        if (!new_data)
            return -4;
        b->data = new_data;
        b->m_data = new_m;
    }
    //if (bgzf_read(fp, b->data, c->l_qname) != c->l_qname) return -4;
    memcpy(b->data, data + cur, c->l_qname);
    cur += c->l_qname;
    for (i = 0; i < c->l_extranul; ++i) b->data[c->l_qname+i] = '\0';
    c->l_qname += c->l_extranul;
    memcpy(b->data + c->l_qname, data + cur, b->l_data - c->l_qname);      
    if (is_be) swap_data(c, b->l_data, b->data, 0);
    return 4 + block_len;
}

ReadBAMChunk* PBBamFile::getOneUncompressedChunk(){
    uint8_t out[CHUNK];
    uint8_t tmp_local[18];  
    ReadBAMChunk* localReadChunk = getOneChunk();
    
    if (bamFile->pbread(tmp_local,1, 18) < 18){        
        return nullptr;
    }
    uint32_t xlen = tmp_local[11] * 256 + tmp_local[10];
    uint32_t bsize = tmp_local[17] * 256 + tmp_local[16];
    int cdata_len = bsize - xlen - 19;        
    localReadChunk->compressedSize = cdata_len;
    assert(localReadChunk->compressedData != nullptr);            
    if (bamFile->pbread(localReadChunk->compressedData, 1, cdata_len) < cdata_len){
        return nullptr;        
    }
    if (bamFile->pbread(tmp_local, 1, 8) < 8){
        return nullptr;
    }

    localReadChunk->offset = 0;
    localReadChunk->uncompressedSize = fromByteArray(tmp_local + 4);       
    assert(localReadChunk->if_processed == false);
    
    uint32_t sizeOut = 0;
    z_stream strm = {0};
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.next_in = localReadChunk->compressedData;
    strm.avail_in = 0;
    CALL_ZLIB (inflateInit2(&strm, -15));

    unsigned have = 0;
    int flush =  Z_NO_FLUSH;
    strm.avail_in = localReadChunk->compressedSize;
    do {                
        have = 0;
        strm.avail_out = CHUNK ;
        strm.next_out = out;              
        CALL_ZLIB (inflate (& strm, flush));                
        have = CHUNK  - strm.avail_out;				
        memcpy(localReadChunk->uncompressedData + sizeOut, out, have);
        sizeOut += have;
    }
    while (strm.avail_out == 0);
    (void)inflateEnd(&strm);
    return localReadChunk;

}
std::string PBBamFile::getString(int len){
    std::string result;
    while(1){
        int bytesNeeded = len - result.length();     
        int bytesAvailable = localSize - localPointer;
        int bytesRead = std::min(bytesAvailable, bytesNeeded);
        if (bytesRead > 0){
            result  += std::string((char *)(localBuffer + localPointer), bytesRead);
            localPointer += bytesRead;
            if (result.length() == len)
                return result;
        }
        
        ReadBAMChunk* tmp = getOneUncompressedChunk();
        assert(tmp != nullptr);
        if (localSize - localPointer > 0)
            memmove(localBuffer, localBuffer + localPointer, localSize - localPointer);
        size_t cur_len = localSize - localPointer;                        
        memcpy(localBuffer + cur_len, tmp->uncompressedData + tmp->offset, tmp->uncompressedSize - tmp->offset);            
        localSize = tmp->uncompressedSize + cur_len - tmp->offset;
        localPointer = 0;
        tmp->if_processed = false;    
        if (ChunkPool.write_available()) {
            assert(ChunkPool.push(tmp));
            chunkPoolSize++;
        }
        else
            delete tmp;            
    }
}
int PBBamFile::getInt(){
    while(1){
        if (localSize > localPointer + 4){
            int result = fromByteArray(localBuffer + localPointer);            
            localPointer += 4;
            return result;

        }
        else{
            ReadBAMChunk* tmp = getOneUncompressedChunk();
            assert(tmp != nullptr);
            if (localSize - localPointer > 0)
                memmove(localBuffer, localBuffer + localPointer, localSize - localPointer);
	        size_t cur_len = localSize - localPointer;                        
            memcpy(localBuffer + cur_len, tmp->uncompressedData + tmp->offset, tmp->uncompressedSize - tmp->offset);            
            localSize = tmp->uncompressedSize + cur_len - tmp->offset;
            localPointer = 0;
            tmp->if_processed = false;    
            if (ChunkPool.write_available()) {
                assert(ChunkPool.push(tmp));
                chunkPoolSize++;
            }
            else
                delete tmp;    
        }
    }
}


// Minimal sanitisation of a header to ensure.
// - null terminated string.
// - all lines start with @ (also implies no blank lines).
//
// Much more could be done, but currently is not, including:
// - checking header types are known (HD, SQ, etc).
// - syntax (eg checking tab separated fields).
// - validating n_targets matches @SQ records.
// - validating target lengths against @SQ records.
 bam_hdr_t *pb_sam_hdr_sanitise(bam_hdr_t *h) {
    if (!h)
        return NULL;

    // Special case for empty headers.
    if (h->l_text == 0)
        return h;

    uint32_t i, lnum = 0;
    char *cp = h->text, last = '\n';
    for (i = 0; i < h->l_text; i++) {
        // NB: l_text excludes terminating nul.  This finds early ones.
        if (cp[i] == 0)
            break;

        // Error on \n[^@], including duplicate newlines
        if (last == '\n') {
            lnum++;
            if (cp[i] != '@') {
                fprintf(stderr, "Malformed SAM header at line %u", lnum);
                bam_hdr_destroy(h);
                return NULL;
            }
        }

        last = cp[i];
    }

    if (i < h->l_text) { // Early nul found.  Complain if not just padding.
        uint32_t j = i;
        while (j < h->l_text && cp[j] == '\0') j++;
        if (j < h->l_text)
            fprintf(stderr,"Unexpected NUL character in header. Possibly truncated");
    }

    // Add trailing newline and/or trailing nul if required.
    if (last != '\n') {
        fprintf(stderr,"Missing trailing newline on SAM header. Possibly truncated");

        if (h->l_text == UINT32_MAX) {
            fprintf(stderr,"No room for extra newline");
            bam_hdr_destroy(h);
            return NULL;
        }

        if (i >= h->l_text - 1) {
            cp = (char *)realloc(h->text, (size_t) h->l_text+2);
            if (!cp) {
                bam_hdr_destroy(h);
                return NULL;
            }
            h->text = cp;
        }
        cp[i++] = '\n';

        // l_text may be larger already due to multiple nul padding
        if (h->l_text < i)
            h->l_text = i;
        cp[h->l_text] = '\0';
    }

    return h;
}
bam_hdr_t * PBBamFile::read_header(){
    bam_hdr_t *h;
    h = bam_hdr_init();
    //if (!h) goto nomem;
    string BAMInitial = getString(4);
    if (BAMInitial.compare(0, 3, "BAM") != 0){
        fprintf(stderr,"Incorrect BAM Header\n");
        exit(EXIT_FAILURE);
    }
    h->l_text = getInt();
    if (is_be) pb_ed_swap_4p(&h->l_text);

    std::string text = getString(h->l_text);
    text += '\0';
    size_t bufsize = ((size_t) h->l_text) + 1;
    if (bufsize < h->l_text) {
        fprintf(stderr,"Large header size error\n");
        exit(EXIT_FAILURE);
    } // so large that adding 1 overflowed
    h->text = (char*)malloc(bufsize);
    memcpy(h->text, text.c_str(), text.length());

    h->n_targets = getInt();
    if (is_be) pb_ed_swap_4p(&h->n_targets);
    if (h->n_targets < 0) {
        fprintf(stderr,"There is no chromosome information in the BAM Header\n");
        exit(EXIT_FAILURE);
    }


    if (h->n_targets > 0) {
        h->target_name = (char**)calloc(h->n_targets, sizeof(char*));
        if (!h->target_name) {
            fprintf(stderr,"Allocation error in the BAM header for target_name\n");
            exit(EXIT_FAILURE);
        } 
        h->target_len = (uint32_t*)calloc(h->n_targets, sizeof(uint32_t));
        if (!h->target_len) {
            fprintf(stderr,"Allocation error in the BAM header for target_len\n");
            exit(EXIT_FAILURE);
        } 
    }
    else {
        h->target_name = NULL;
        h->target_len = NULL;
    }
    for (int i = 0; i < h->n_targets; i++){
        
        int l_name = getInt();
        if (is_be) pb_ed_swap_4p(&l_name);
        if (l_name <= 0)  {
            fprintf(stderr,"Name length of chromosome number %d in the BAM header is negative\n", l_name);
            exit(EXIT_FAILURE);
        } 
        std::string name = getString(l_name);
        if (name[name.length() -1] != '\0'){
            name += '\0';
        }        
        h->target_name[i] = (char*)malloc(name.length());        
        memcpy(h->target_name[i], name.c_str(), name.length());
        h->target_len[i] = getInt();
        if (is_be) pb_ed_swap_4p(&h->target_len[i]);
    }
   
    return pb_sam_hdr_sanitise(h);
   
}

bool PBBamFile::read_bam(bam1_t * b){    
    assert(b!=nullptr);
    while(1){
        int ret = -1;
        if (localSize > localPointer){
            ret  = bam_read1_pb(b, localBuffer + localPointer, localSize - localPointer, false);
        }
        if (ret > 0){	            
            auto pos = b->core.pos + 1;
            auto end = bam_endpos(b) ;
            auto tid = b->core.tid;
            localPointer += ret; 
       
            assert(b!=nullptr);
            return true;
                          
        } else{
            if (localSize - localPointer > 0)
                memmove(localBuffer, localBuffer + localPointer, localSize - localPointer);
	        size_t cur_len = localSize - localPointer;            
            ReadBAMChunk* tmp = nullptr;
            while ((!read_done_local)  && (!ChunkReadQueue.pop(tmp)));
            if (tmp == nullptr){
                if (!ChunkReadQueue.pop(tmp)){

                    return false;
                }
            }

            assert(tmp != nullptr);

            while(tmp->if_processed == false);
            localSize = tmp->uncompressedSize + cur_len - tmp->offset;
            if (localSize > localBufferCapacity){
                while(localSize > localBufferCapacity)
                    localBufferCapacity *= 2;
                localBuffer = (uint8_t*)realloc(localBuffer, localBufferCapacity);
            }
            memcpy(localBuffer + cur_len, tmp->uncompressedData + tmp->offset, tmp->uncompressedSize - tmp->offset);
            
            localPointer = 0;
            tmp->if_processed = false;
            if (ChunkPool.write_available()) {
                assert(ChunkPool.push(tmp));
                chunkPoolSize++;
            }
            else
                delete tmp;
        }
    }

}

void getGZReadThread_parallel(ChunkLockFreeQueue& ChunkUnzipQueue, std::atomic<uint32_t>& UnzipQueueSize, std::atomic<bool>& read_done_local){

	ReadBAMChunk* tmp;
	uint8_t* dataIn_local;
	uint8_t* dataOut_local;
	uint32_t dataSize;
	uint8_t out[CHUNK];

	while((!read_done_local) || (UnzipQueueSize > 0) ){
		if (!ChunkUnzipQueue.pop(tmp))
			continue;
		--UnzipQueueSize;

		uint32_t sizeOut = 0;
		z_stream strm = {0};
		strm.zalloc = Z_NULL;
		strm.zfree = Z_NULL;
		strm.opaque = Z_NULL;
		strm.next_in = tmp->compressedData;
		strm.avail_in = 0;
		CALL_ZLIB (inflateInit2(&strm, -15));

		unsigned have = 0;
		int flush =  Z_NO_FLUSH;
		strm.avail_in = tmp->compressedSize;
		do {                
			have = 0;
			strm.avail_out = CHUNK ;
			strm.next_out = out;              
			CALL_ZLIB (inflate (& strm, flush));                
			have = CHUNK  - strm.avail_out;				
			memcpy(tmp->uncompressedData + sizeOut, out, have);
			sizeOut += have;
		}
		while (strm.avail_out == 0);
		(void)inflateEnd(&strm);

		tmp->if_processed = true;
	}         
}

ReadBAMChunk* PBBamFile::getOneChunk(){
    if (chunkPoolSize > 0) {
		ReadBAMChunk* tmp;
		assert(ChunkPool.pop(tmp));
        chunkPoolSize--;
        assert(tmp->compressedData != nullptr);
		return tmp;
	}
	if (backupIdx == STEP_MALLOC_SIZE){
		backupData = (uint8_t*) malloc(STEP_MALLOC_SIZE * BLOCK_SIZE * sizeof(uint8_t));
		backupUncompressedData = (uint8_t*) malloc(STEP_MALLOC_SIZE * BLOCK_SIZE * sizeof(uint8_t));
		freeQueue.push_back(backupData);
		freeQueue.push_back(backupUncompressedData);
		backupIdx = 0;
	}
	ReadBAMChunk* tmp = new ReadBAMChunk(backupIdx * BLOCK_SIZE * sizeof(uint8_t) + backupData, backupIdx * BLOCK_SIZE * sizeof(uint8_t) + backupUncompressedData);
    assert(tmp->compressedData != nullptr);
	++backupIdx;
	return tmp;
}

void PBBamFile::fillChunks_stream(std::atomic<bool>& read_done_local, PBLocalFile * bamFile){
    uint8_t tmp_local[18];
    PBLoopReader * pbloopReader = new PBLoopReader(bamFile);

    uint64_t cur_chunk = 0;
    pbloopReader->reset();       
    size_t index = 0;
    while(!read_done_local){			
        if (pbloopReader->copy(tmp_local,18) < 18){
            
            break;
        }
        uint32_t xlen = tmp_local[11] * 256 + tmp_local[10];
        uint32_t bsize = tmp_local[17] * 256 + tmp_local[16];

        int cdata_len = bsize - xlen - 19;

        ReadBAMChunk* localReadChunk = getOneChunk();
        localReadChunk->compressedSize = cdata_len;
        assert(localReadChunk->compressedData != nullptr);            
        pbloopReader->copy(localReadChunk->compressedData, cdata_len);        
        pbloopReader->copy(tmp_local, 8);
        
        localReadChunk->offset = 0;
        localReadChunk->uncompressedSize = fromByteArray(tmp_local + 4);       
        assert(localReadChunk->if_processed == false);
        while(ChunkReadQueue.write_available() == 0);
        assert(ChunkReadQueue.push(localReadChunk));
        while(ChunkUnzipQueue.push(localReadChunk) == false);
        ++UnzipQueueSize;
        cur_chunk += (bsize + 1);

    }    
            
   
    delete pbloopReader;
    read_done_local = true;
}

PBBamFile::PBBamFile(const char * _filename, int _nZipThreads, const char * mode, uint32_t padding, int intervalPadding, bool longRead):
                    backupData(nullptr), backupUncompressedData(nullptr), backupIdx(STEP_MALLOC_SIZE),
					ChunkUnzipQueue(longRead ? POOL_SIZE_LONG_READ : POOL_SIZE),  UnzipQueueSize(0),
                    ChunkPool(longRead ? POOL_SIZE_LONG_READ : POOL_SIZE), ChunkReadQueue(longRead ? POOL_SIZE_LONG_READ : POOL_SIZE) {
    filename = _filename;
    readingMode = false;

    chunkPoolSize = 0;
    assert(mode[0] == 'r');
    readingMode = true;
    is_be = ed_is_big();  
    bamFile = new PBLocalFile(filename, mode);
    if (bamFile->isFailed()){
        cerr << "Cannot open file:" << filename << ". Exiting..." << endl;
        exit(EXIT_FAILURE);
    }

    localBuffer = new uint8_t[2 * BLOCK_SIZE];
    localBufferCapacity = 2 * BLOCK_SIZE;
    localSize = 0;
    localPointer = 0;


    nZipThreads = _nZipThreads;

    uint8_t* data = (uint8_t*) malloc(INIT_CHUNK_QUEUE_SIZE * BLOCK_SIZE * sizeof(uint8_t));
    uint8_t* output = (uint8_t *) malloc(INIT_CHUNK_QUEUE_SIZE * BLOCK_SIZE * sizeof(uint8_t));

    for (uint32_t i = 0; i < INIT_CHUNK_QUEUE_SIZE; ++i){
        assert(ChunkPool.push(new ReadBAMChunk(i * BLOCK_SIZE * sizeof(uint8_t) + data, i * BLOCK_SIZE * sizeof(uint8_t) + output)));
        chunkPoolSize++;
    }
    freeQueue.push_back(data);
    freeQueue.push_back(output);
    mHeader  = read_header();
    read_done_local = false;


    fillBamThread = PBThread("fillbamStream", &PBBamFile::fillChunks_stream, this ,std::ref(read_done_local), bamFile);

    unzip_threads = new PBThread[nZipThreads];
    for (int i = 0; i < nZipThreads; ++i)
        unzip_threads[i] = PBThread("UnZipBam_" + std::to_string(i), getGZReadThread_parallel, std::ref(ChunkUnzipQueue), std::ref(UnzipQueueSize), std::ref(read_done_local));

    
     
}

PBBamFile::~PBBamFile(){
   assert(readingMode);
    read_done_local = true;     
    for (int i = 0; i < nZipThreads; ++i)
	    unzip_threads[i].join();
    ReadBAMChunk* readChunk;        
    while (ChunkReadQueue.pop(readChunk)) 
        delete readChunk;
    fillBamThread.join();
    delete [] localBuffer;
    delete bamFile;
    for (auto x : freeQueue)
        free(x);
    while (ChunkReadQueue.pop(readChunk)) 
        delete readChunk;
    
    while (ChunkUnzipQueue.pop(readChunk)) {
        --UnzipQueueSize;
    }
    while (ChunkPool.pop(readChunk)){
        chunkPoolSize--;
        delete readChunk;
    }

    bam_hdr_destroy(mHeader);
    delete [] unzip_threads;
    
}

#define LOOPSIZE 10000000
PBLoopReader::PBLoopReader(PBLocalFile * _fp){
    fp = _fp;
    position = 0;
    size = 0;
    data = new uint8_t[LOOPSIZE];
}
void PBLoopReader::reset(){
    position = 0;
    size = 0;
}

PBLoopReader::~PBLoopReader(){
    delete [] data;
}
int PBLoopReader::copy(uint8_t * ptr, uint count){    
    assert (count < LOOPSIZE);
    if (size - position >= count){
        memcpy(ptr, data + position, count);
        position += count;
        return count;
    }
    else{
        if (size != position)
            memmove(data, data + position, size - position);
        int readSize = fp->pbread(data + size - position, 1,  LOOPSIZE - size + position);
        if ( readSize < 0){
            if (size <= position)
                return -1;
            memcpy(ptr, data, size - position);
            int readSize = size - position;
            position = size - position;
            size = position;
            return readSize;
        }else{
            size = readSize + size - position;
            position = 0; 
            readSize = std::min(count, size);
            memcpy(ptr, data + position, readSize);
            position += readSize;
            return readSize;
        }
    }
}
