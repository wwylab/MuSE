

#ifndef __PBBAMFILE_H__
#define __PBBAMFILE_H__
#include <stdio.h>

#include <boost/lockfree/queue.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#include <vector>
#include <queue>

#include "PBLocalFile.h"
#include "pbthreads.h"

#include "htslib/sam.h"


#define CHUNK 100000
#define POOL_SIZE 100000
#define POOL_SIZE_LONG_READ 2000
#define INIT_CHUNK_QUEUE_SIZE 2000
#define STEP_MALLOC_SIZE 100
#define BLOCK_SIZE 0x10000

struct ReadBAMChunk{
public:
	uint8_t* compressedData;
	uint8_t* uncompressedData;
	uint16_t compressedSize;
	uint32_t uncompressedSize;
	uint32_t offset;
	std::atomic<bool> if_processed;

	ReadBAMChunk(uint8_t* compressed_in, uint8_t* uncompressed_in):compressedData(compressed_in), uncompressedData(uncompressed_in),
																compressedSize(0), uncompressedSize(0), if_processed(false){};
};

class PBLoopReader
{
    private:
        uint8_t * data;
        PBLocalFile * fp;
        uint32_t position;
        uint32_t size;
    public:
        void reset();
        int copy(uint8_t * ptr, uint count);
        PBLoopReader(PBLocalFile * _fp);
        ~PBLoopReader();
};

typedef boost::lockfree::queue<ReadBAMChunk*> ChunkLockFreeQueue;
typedef boost::lockfree::spsc_queue<ReadBAMChunk*> ChunkSpscQueue;


class PBBamFile
{
    private:
        bool readingMode;
        int nZipThreads;
        const char * filename;

        std::atomic_bool readingDone;
        std::atomic<bool> read_done_local;

        bool is_be;

        uint8_t* localBuffer;
        uint64_t localBufferCapacity;
        uint32_t localPointer;
        uint32_t localSize;

       	ChunkLockFreeQueue ChunkUnzipQueue;
        std::atomic<uint32_t> UnzipQueueSize;
        int backupIdx;	
        uint8_t* backupData;
	    uint8_t* backupUncompressedData;
	    std::vector<uint8_t*> freeQueue;

        PBLocalFile * bamFile;

        PBThread*  unzip_threads;
        PBThread fillBamThread;
               
        ChunkSpscQueue	ChunkPool;
        std::atomic_int chunkPoolSize;
	    ChunkSpscQueue	ChunkReadQueue;

        ReadBAMChunk* getOneChunk();
        void fillChunks_stream(std::atomic<bool>& read_done_local, PBLocalFile * bamFile);

    public:
        PBBamFile(const char * _filename, int _nZipThreads, const char * mode, uint32_t padding = 100, int intervalPadding = 0, bool longRead = false);
        ~PBBamFile(); 
        bam_hdr_t *  mHeader;   

        bool read_bam(bam1_t*);   
        bam_hdr_t * read_header();     
        std::string getString(int len);
        int getInt();
        ReadBAMChunk* getOneUncompressedChunk();

        inline int getTid(const char *ref){return bam_name2id(mHeader, ref); }
       	inline uint32_t getChrLength(int32_t chrTid) { return mHeader->target_len[chrTid]; }
        bam_hdr_t* getHeader(){return mHeader;}
};



#endif
