#ifndef __MUSE_READER_H__
#define __MUSE_READER_H__

#include <boost/lockfree/queue.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#include <list>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <queue>
#include "muse_const.h"
#include "PBBamFile.h"

extern std::string Version;
extern std::string buildDate;  // e.g. 'Dec 15 2009'
extern std::string buildTime;  // e.g. '15:25:56'


class Region;
#define PILEUP_SIZE 10000
#define REGION_SIZE 1000

using PileUpLockFreeQ = boost::lockfree::queue<Region*>;
using PileupSpscQ = boost::lockfree::spsc_queue<Region*, boost::lockfree::capacity<PILEUP_SIZE>>;

using READQUEUETYPE = boost::lockfree::queue<bam1_t*>;
using READQUEUESPSCTYPE = boost::lockfree::spsc_queue<bam1_t*, boost::lockfree::capacity<200000> >;


struct node{
	bam1_t* b;
	int64_t end;
};

class PBReader {
public:
	PBBamFile* rFP_tumor, *rFP_normal;
	bam_hdr_t* mHeader_tumor, *mHeader_normal;
	std::atomic<bool> ReadingDone;
	READQUEUESPSCTYPE ReadQueue;
	READQUEUESPSCTYPE ReadPool;
	bam1_t* curPtr = nullptr;
	std::multimap<uint64_t, bam1_t*> mergedReads;
	mplp_conf_t* conf;

	const char* chrName(int i){
		if (i >= mHeader_normal->n_targets){
			std::cerr << "Wrong chromosome idx " << i << std::endl;
			exit(EXIT_FAILURE);
		}
		return mHeader_normal->target_name[i];
	}
	PBReader(std::string& tumorName, std::string& normalName, mplp_conf_t* conf_in): conf(conf_in){
		ReadingDone.store(false);
		curPtr = nullptr;   	

		for (auto& fileName : {tumorName, normalName}){
			auto pos = fileName.find_last_of('.');
			if (pos == std::string::npos || fileName.substr(pos) != ".bam"){
				std::cerr << "Input file should have extension .bam" << std::endl;
				exit(EXIT_FAILURE);
			}
		}
		rFP_tumor = new PBBamFile(tumorName.c_str(), 2, "r", 0);
		mHeader_tumor = rFP_tumor->getHeader();
		rFP_normal = new PBBamFile(normalName.c_str(), 2, "r", 0);
		mHeader_normal = rFP_normal->getHeader();

		
		assert(!(conf->flag & MPLP_ILLUMINA13));
	}

	~PBReader(){
		delete rFP_normal; delete rFP_tumor;
		if (curPtr != nullptr)            	
			bam_destroy1(curPtr);
		while(ReadPool.read_available()!= 0){
			bam1_t * read;
			read = ReadPool.front();
			ReadPool.pop();
			bam_destroy1(read);
		}
	}

	void setDone(){ ReadingDone.store(true);}
	int getReadQueueSize() { return ReadQueue.read_available(); }
	size_t getReadPoolSize() { return ReadPool.read_available();}

	void enqueue(bam1_t * bamPtr){
		assert(bamPtr);
		while(ReadQueue.write_available() == 0);
		assert(ReadQueue.push(bamPtr));
	}
	
	bam1_t* getNextFilteredRecordFromQueue(){
		bam1_t* result;
		while(!ReadingDone){
			if (ReadQueue.pop(result)){
				return result;
			}
		}

		if (ReadQueue.pop(result)){
			return result;
		} else {
			return nullptr;
		}	
	}
	
	bam1_t* getBam(){
		bam1_t* curPtr;
		if (ReadPool.read_available() == 0)
			curPtr = bam_init1();
		else{
			curPtr = ReadPool.front();
			ReadPool.pop();
		}
		return curPtr;
	}
	
	void release(bam1_t* read){
		if (ReadPool.write_available() == 0)
			bam_destroy1(read);
		else	
			ReadPool.push(read);
	}

	bool applyFilters(bam1_t* b){
		if (b->core.flag&BAM_FUNMAP)
			return false;
		if (b->core.flag & conf->flag_mask)
			return false;
		return true;
	}

	void fillBamQueue(){
		curPtr = getBam();
		std::vector<PBBamFile*> rFP = {rFP_tumor, rFP_normal};
		for (int i = 0; i < rFP.size(); ++i){
			auto& x = rFP[i];
			while(x->read_bam(curPtr)){		
				if (curPtr->core.tid == -1){
					break;
				}
				if(applyFilters(curPtr)) {
					curPtr->core.unused1 = i;
					mergedReads.insert(std::make_pair(uint64_t(curPtr->core.tid) << 32 | curPtr->core.pos, curPtr));
					curPtr = getBam();
					break;
				}		                 
			}
		}

		while(mergedReads.size()){
			auto read = mergedReads.begin()->second;
			mergedReads.erase(mergedReads.begin());
			uint8_t fileIdx = read->core.unused1;
			enqueue(read);
			while(rFP[fileIdx]->read_bam(curPtr)){
				if (curPtr->core.tid == -1){
					break;
				}
				if (applyFilters(curPtr)){
					curPtr->core.unused1 = fileIdx;
					mergedReads.insert(std::make_pair(uint64_t(curPtr->core.tid) << 32 | curPtr->core.pos, curPtr));
					curPtr = getBam();
					break;
				}
			}
		}

	}

};


class Region{
public:
	std::vector<node> nodes;
	int tid;
	int64_t start, end;
	std::atomic<bool> ifDone;
	std::string output;

	Region():ifDone(false){}
	void init(std::deque<node>& overlapQueue, const node& Node, int& lastRegionTid, int64_t& lastRegionPos){
		
		bam1_t* b;
		if (overlapQueue.size()){
			b = overlapQueue.front().b;
		}else{
			b = Node.b;
		}
		tid = b->core.tid;
		start = (b->core.pos / REGION_SIZE) * REGION_SIZE;
		if (tid == lastRegionTid && start <= lastRegionPos){
			start = lastRegionPos + REGION_SIZE;
		}
		end = start + REGION_SIZE - 1;
		lastRegionTid = tid;
		lastRegionPos = start;
		
		while(overlapQueue.size()){
			auto& top = overlapQueue.front();
			
			if (top.b->core.pos == top.end){
				overlapQueue.pop_front();
				continue;
			}

			if (tid == top.b->core.tid && top.b->core.pos <= end && start <= top.end - 1){
				nodes.push_back(top);
			}else{
				break;
			}
			overlapQueue.pop_front();
		}
		if (Node.b){
			if(Node.b->core.tid == tid && Node.b->core.pos <= end && start <= Node.end - 1){
				nodes.push_back(Node);
			}else{
				overlapQueue.push_back(Node);
			}
		}
	}
	void copyBackToQueue(std::deque<node>& overlapQueue){
		for (int i = nodes.size() - 1; i >= 0; --i){
			if (nodes[i].end > end + 1){
				overlapQueue.push_front(nodes[i]);
			}
		}

	}
};

struct nodeCmp{
	bool operator()(const std::pair<std::pair<int, int64_t>, std::list<std::pair<int, int>>::iterator>& left, 
					const std::pair<std::pair<int, int64_t>, std::list<std::pair<int, int>>::iterator>& right) const{
		if (left.first.first != right.first.first)
			return left.first.first > right.first.first;
		return left.first.second > right.first.second;
	}
};

class bam_plp{
public:
	std::list<std::pair<int, int>> pendingReadsPosAll[2];
	std::priority_queue<std::pair<std::pair<int, int64_t>, std::list<std::pair<int, int>>::iterator>, 
						std::vector<std::pair<std::pair<int, int64_t>, std::list<std::pair<int, int>>::iterator>>,
						nodeCmp> pendingReads[2];
	std::deque<node> overlapQueue;
	int lastRegionTid;
	int64_t lastRegionPos;
	int tid[2];
	int64_t pos[2];
	Region* region;
	int max_cnt;

	bam_plp(int max_cnt_in):lastRegionTid(-1), lastRegionPos(-1), region(nullptr), max_cnt(max_cnt_in){
		tid[0] = tid[1] = -1;
		pos[0] = pos[1] = -1;
	};
	Region* push(node& Node){

		auto& tid_local = tid[Node.b->core.unused1];
		auto& pos_local = pos[Node.b->core.unused1];
		auto& pendingRead = pendingReads[Node.b->core.unused1];

		auto& pendingReadsPos = pendingReadsPosAll[Node.b->core.unused1];

		if (pendingRead.size() + 2 > max_cnt && Node.b->core.tid == tid_local && Node.b->core.pos == pos_local){
			bam_destroy1(Node.b);
			return nullptr;
		}

		int max_tid = Node.b->core.tid;
		int64_t max_pos = Node.b->core.pos;

		pendingReadsPos.emplace_back(Node.b->core.tid, Node.b->core.pos);
		pendingRead.push(std::make_pair(std::make_pair(Node.b->core.tid, Node.end), --pendingReadsPos.end()));

		Region* toReturn = nullptr;
		if (region){
			if (Node.b->core.tid > region->tid || (Node.b->core.tid == region->tid && Node.b->core.pos > region->end + 1)){
				region->copyBackToQueue(overlapQueue);
				toReturn = region;
				region = nullptr;
			}
		}

		if (!region){
			region = new Region();
			region->init(overlapQueue, Node, lastRegionTid, lastRegionPos);
		}else{
			
			if (Node.b->core.pos <= region->end)
				region->nodes.push_back(Node);
			else{
				overlapQueue.push_back(Node);
			}
		}

		while(pendingRead.size()){
			auto& tmpNode = pendingRead.top();
			if (tmpNode.first.first < Node.b->core.tid || tmpNode.first.second <= Node.b->core.pos - 1){
				pendingReadsPos.erase(tmpNode.second);
				pendingRead.pop();
			}else{
				break;
			}
		}

		
		while(tid_local < max_tid || (tid_local == max_tid && pos_local < max_pos)){
			assert(pendingReadsPos.size());
			int tmp_tid = pendingReadsPos.front().first;
			int tmp_pos = pendingReadsPos.front().second;
			if (tid_local < tmp_tid){
				tid_local = tmp_tid;
				pos_local = tmp_pos;
			}else if(pos_local < tmp_pos){
				pos_local = tmp_pos;
			}else
				++pos_local;
		}

		return toReturn;
	}

	Region* pop(){
		Region* toReturn = nullptr;
		if (region){
			region->copyBackToQueue(overlapQueue);
			toReturn = region;
			region = nullptr;
		}

		if (overlapQueue.size()){
			region = new Region();
			region->init(overlapQueue, {nullptr, 0}, lastRegionTid, lastRegionPos);
		}
		return toReturn;
	}

};




void construct_convert_table();
void WriteHeader(PBLocalFile* outFile, int argc, char * const argv[], bam_hdr_t *tumorH, bam_hdr_t *normalH, const std::string& refGenome);
void processPileup(PBReader* reader, mplp_conf_t* conf, std::atomic<bool>& gatherDone, PileUpLockFreeQ& processQ, std::atomic<uint32_t>& processQSize);

#endif