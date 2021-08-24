

#ifndef __MUSE_WORKERMEM_H__
#define __MUSE_WORKERMEM_H__

#include "muse_const.h"
#include <unordered_set>
#include <queue>
#include <list>

//charArrayHashFun: cited from http://www.cse.yorku.ca/~oz/hash.html
size_t charArrayHashFun(const char* s){
	size_t h = 5381;
	int c;
	while ((c = *s++))
		h = ((h << 5) + h) + c;
	return h;
}

struct charArrayHash{
	size_t operator()(const std::pair<const char*, size_t>& item) const{
		return item.second;
	}
};

struct charArrayEuqal{
	bool operator()(const std::pair<const char*, size_t>& left, const std::pair<const char*, size_t>& right) const{
		return (strcmp(left.first, right.first) == 0);
	}
};

struct cmp{
	bool operator()(const std::pair<int, std::list<bam_pileup1_t_pb*>::iterator>& left, 
					const std::pair<int, std::list<bam_pileup1_t_pb*>::iterator>& right) const{
		return left.first > right.first;
	}
};

class WorkMem{
public:
	std::ostringstream os;
	std::vector<bam_pileup1_t_pb*> plpToFree;
	std::vector<bam_pileup1_t_pb*> all_plp;
	
	std::vector<std::vector<bam_pileup1_t_pb*>> plp;
	std::vector<std::string> rawS, s;
	std::vector<std::vector<int>> rawQ, q;
	std::unordered_map<std::pair<const char*, size_t>, std::pair<int, int>, charArrayHash, charArrayEuqal> overlapPair[2];
	std::unordered_set<std::string> overlapRG;
	std::vector<bam_pileup1_t_pb*> rawTumorPileupElementPointer;
	std::vector<bam_pileup1_t_pb*> noOverlapNoDelPassMinQNoMapQ0TumorPileupElementPointer;
	std::vector<bam_pileup1_t_pb*> finalNormalPileupElementPointer;

	std::pair<int, int> tumorVariantAlleleCount[4];
	std::pair<int, int> normalVariantAlleleCountAndQuality[4];
	std::pair<int,int> rawNormalAlleleCountVec[5];
	int rawNormalAlleleCount[5];
	int rawNormalAlleleQuality[5];
	int rawTumorAlleleCount[5];
	int rawTumorAlleleQuality[5];

	std::vector<double> smallestDistance;

	std::map<int, int> QualityCount[4];

	std::vector<double> unscaledBaseFreq;
	double brlensMin;
	double brlensMax;
	double min_diff_lk_local;

	std::vector<int>    weight;
	std::vector<double> condLike;
	std::vector<double> pMatrix;
	std::vector<double> baseFreq;	// if BAYESIAN, prior is Dirichlet(1,1,1,1)
	double              kappa;		// if BAYESIAN, prior is Exp(1)
	double              brlens;		// if BAYESIAN, prior is Exp(brlenspr_mean)
	double              lnlike;
	char                refLetter;
	char                altLetter;
	bool   optimizeBaseFreq;

	std::string altString;
	std::string normalFormat;
	std::string tumorFormat;

	std::vector<std::vector<int>> plp_push_idx;
	std::vector<std::vector<int>> plp_pop_idx;

	WorkMem():plp(2), rawS(2), rawQ(2), s(2), q(2), plp_push_idx(REGION_SIZE), plp_pop_idx(REGION_SIZE){
		os.precision(5);
		os << std::scientific;
		bam_pileup1_t_pb* plpMem = (bam_pileup1_t_pb*)malloc(10000 * sizeof(bam_pileup1_t_pb));
		plpToFree.push_back(plpMem);
		for (int i = 0; i < 10000; ++i)
			all_plp.push_back(plpMem + i);
	}

	~WorkMem(){
		for (auto x: plpToFree)
			free(x);
	}
	void clear(){
		for (auto& x: plp)
			x.clear();
	}

	std::list<bam_pileup1_t_pb*> plp_list;
	std::priority_queue<std::pair<int, std::list<bam_pileup1_t_pb*>::iterator>, std::vector<std::pair<int, std::list<bam_pileup1_t_pb*>::iterator>>, cmp> plp_q;
	void pareparePlp(std::vector<node>& nodes){
		if (nodes.size() > all_plp.size()){
			auto allocNum = nodes.size() - all_plp.size();
			bam_pileup1_t_pb* tmp = (bam_pileup1_t_pb*)malloc(allocNum * sizeof(bam_pileup1_t_pb));
			plpToFree.push_back(tmp);
			for (int j = 0; j < allocNum; ++j)
				all_plp.push_back(tmp + j);
		}
		
		plp_list.clear();
		decltype(plp_q)().swap(plp_q);
		for (int i = 0; i < nodes.size(); ++i){
			all_plp[i]->b = nodes[i].b;
			all_plp[i]->end = nodes[i].end;
			all_plp[i]->s.k = 0;
			plp_list.push_back(all_plp[i]);
			plp_q.push(make_pair(nodes[i].end, --plp_list.end()));
		}
	}
	bool getPlp(int pos){
		while(plp_q.size() && plp_q.top().first <= pos){
			plp_list.erase(plp_q.top().second);
			plp_q.pop();
		}

		for (auto iter = plp_list.begin(); iter != plp_list.end(); ++iter){

			bam_pileup1_t_pb* tmp = *iter;
			if (tmp->b->core.pos > pos)
				break;
			
			plp[tmp->b->core.unused1].push_back(tmp);
		}
		return plp[0].size() && plp[1].size();
	}
	void processPileupFun(bam_hdr_t* h, mplp_conf_t* conf, std::string& output, int tid, uint32_t pos, const char* ref, uint32_t ref_len);
	void ReformatData(std::string &s, std::vector<int> &q, std::vector<int> &weight, std::vector<double> &condLike);
	void ResetParams();
	void Check_Br_Len_Bounds();
	void NormalizeBaseFreq();
	void CalcProbMat(std::vector<double> &pMat, const std::vector<double> &state_freqs, const double &kappa, const double &brlens);
	void CalcLnLike(double &lnlike, const std::vector<int> &weight, const std::vector<double> &condLike, const std::vector<double> &baseFreq, const std::vector<double> &pMatrix);
	void Adjust_Min_Diff_Lk();
	double Lk();
	double Generic_Brent_Lk(double *param, double ax, double cx, double tol, int n_iter_max, int quickdirty);
	double Br_Len_Brent(double *brlens);
	void Optimize_Br_Len_Serie();
	void Optimize_State_Freqs();
	void Optimiz_All_Free_Param();
	void Round_Optimize(int n_round_max);
	void CalcEvolDistance(std::string &s, std::vector<int> &q);
	std::string scientificStr(double num);
};


#endif