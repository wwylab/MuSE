#include "muse_reader.h"
#include "muse_workerMem.h"
#include "htslib/kstring.h"
#include <cmath>
#include <cfloat>

using namespace std;

std::string Version = "v2.0";
std::string buildDate = __DATE__;  // e.g. 'Dec 15 2009'
std::string buildTime = __TIME__;  // e.g. '15:25:56'


const char *bam_nt16_rev_table = "=ACMGRSVTWYHKDBN";

unsigned char bam_nt16_table[256] = {
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	 1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

char bam_nt16_nt4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };

uint8_t convert_table[256];
char revert_table[5];

void construct_convert_table(){
	memset(convert_table, 0, 256);
	
	convert_table['C'] = 1;
	convert_table['G'] = 2;
	convert_table['T'] = 3;
	convert_table['?'] = 4;

	revert_table[0] = 'A';
	revert_table[1] = 'C';
	revert_table[2] = 'G';
	revert_table[3] = 'T';
	revert_table[4] = '?';
}

inline void resolve_cigar2_pb2(bam_pileup1_t_pb *p, uint32_t pos, cstate_t& s){


	bam1_t *b = p->b;
	bam1_core_t *c = &b->core;
	uint32_t *cigar = bam_get_cigar(b);

	if (s.k == 0){
		s.x = c->pos;
		s.y = 0;
	}
	for (; s.k < c->n_cigar; ++s.k){
		int op = _cop(cigar[s.k]);
		int l = _cln(cigar[s.k]);
		if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP ||
			op == BAM_CEQUAL || op == BAM_CDIFF){
			if (s.x + l > pos)
				break;
			s.x += l;
			if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)
				s.y += l;
		}else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
			s.y += l;
		}
	}

	// collect pileup information
	int op, l;
	op = _cop(cigar[s.k]); l = _cln(cigar[s.k]);
	p->is_del = p->indel = p->is_refskip = 0;
	if (s.x + l - 1 == pos && s.k + 1 < c->n_cigar) { // peek the next operation
		int op2 = _cop(cigar[s.k+1]);
		int l2 = _cln(cigar[s.k+1]);
		if (op2 == BAM_CDEL) p->indel = -(int)l2;
		else if (op2 == BAM_CINS) p->indel = l2;
		else if (op2 == BAM_CPAD && s.k + 2 < c->n_cigar) { // no working for adjacent padding
			int l3 = 0;
			for (int k = s.k + 2; k < c->n_cigar; ++k) {
				op2 = _cop(cigar[k]); l2 = _cln(cigar[k]);
				if (op2 == BAM_CINS) l3 += l2;
				else if (op2 == BAM_CDEL || op2 == BAM_CMATCH || op2 == BAM_CREF_SKIP || op2 == BAM_CEQUAL || op2 == BAM_CDIFF) break;
			}
			if (l3 > 0) p->indel = l3;
		}
	}
	if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
		p->qpos = s.y + (pos - s.x);
	} else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
		p->is_del = 1; 
		p->qpos = s.y; // FIXME: distinguish D and N!!!!!
		p->is_refskip = (op == BAM_CREF_SKIP);
	} // cannot be other operations; otherwise a bug
	p->is_head = (pos == c->pos); p->is_tail = (pos == p->end - 1);
	
	//p->cigar_ind = s->k;
}

bool pileup_seq(const bam_pileup1_t_pb *p, int rb, std::string &s) {
	if(!p->is_del) {
		int c = bam_nt16_rev_table[bam_seqi(bam_get_seq(p->b), p->qpos)];
		//int rb = pos < ref_len ? (*ref_ptr)[pos] : 'N';
		if(c=='=' || bam_nt16_table[c]==bam_nt16_table[rb]) {
			//if(rb == 'N')
			//	s += '?';
			//else
			s.push_back(char((rb)));
			return false;
		}
		else {
			if(c == 'N'){
				s.push_back('?');
				return false;
			}
			else{
				s.push_back(char((c)));
				return true;
			}
		}
	}
	else {
		s.push_back('?');
		return false;
	}
}

// log\binom{n}{k}
inline double lbinom(int n, int k)
{
	if (k == 0 || n == k) return 0;
	return lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
}

// n11  n12  | n1_
// n21  n22  | n2_
//-----------+----
// n_1  n_2  | n

// hypergeometric distribution
inline double hypergeo(int n11, int n1_, int n_1, int n)
{
	return exp(lbinom(n1_, n11) + lbinom(n-n1_, n_1-n11) - lbinom(n, n_1));
}

struct hgacc_t{
	int n11, n1_, n_1, n;
	double p;
};

// incremental version of hypergenometric distribution
double hypergeo_acc(int n11, int n1_, int n_1, int n, hgacc_t *aux)
{
	if (n1_ || n_1 || n) {
		aux->n11 = n11; aux->n1_ = n1_; aux->n_1 = n_1; aux->n = n;
	} else { // then only n11 changed; the rest fixed
		if (n11%11 && n11 + aux->n - aux->n1_ - aux->n_1) {
			if (n11 == aux->n11 + 1) { // incremental
				aux->p *= (double)(aux->n1_ - aux->n11) / n11
					* (aux->n_1 - aux->n11) / (n11 + aux->n - aux->n1_ - aux->n_1);
				aux->n11 = n11;
				return aux->p;
			}
			if (n11 == aux->n11 - 1) { // incremental
				aux->p *= (double)aux->n11 / (aux->n1_ - n11)
					* (aux->n11 + aux->n - aux->n1_ - aux->n_1) / (aux->n_1 - n11);
				aux->n11 = n11;
				return aux->p;
			}
		}
		aux->n11 = n11;
	}
	aux->p = hypergeo(aux->n11, aux->n1_, aux->n_1, aux->n);
	return aux->p;
}

void kt_fisher_exact(int n11, int n12, int n21, int n22, double *two){
	int i, j, max, min;
	double p, q, left, right;
	hgacc_t aux;
	int n1_, n_1, n;

	n1_ = n11 + n12; n_1 = n11 + n21; n = n11 + n12 + n21 + n22; // calculate n1_, n_1 and n
	max = (n_1 < n1_) ? n_1 : n1_; // max n11, for right tail
	min = n1_ + n_1 - n;
	if (min < 0) min = 0; // min n11, for left tail
	//*two = *_left = *_right = 1.;
	*two = 1.;
	if (min == max) return; // no need to do test
	q = hypergeo_acc(n11, n1_, n_1, n, &aux); // the probability of the current table
	// left tail
	p = hypergeo_acc(min, 0, 0, 0, &aux);
	for (left = 0., i = min + 1; p < 0.99999999 * q; ++i) // loop until underflow
		left += p, p = hypergeo_acc(i, 0, 0, 0, &aux);
	--i;
	if (p < 1.00000001 * q) left += p;
	else --i;
	// right tail
	p = hypergeo_acc(max, 0, 0, 0, &aux);
	for (right = 0., j = max - 1; p < 0.99999999 * q; --j) // loop until underflow
		right += p, p = hypergeo_acc(j, 0, 0, 0, &aux);
	++j;
	if (p < 1.00000001 * q) right += p;
	else ++j;
	// two-tail
	*two = left + right;
	if (*two > 1.) *two = 1.;
	// adjust left and right
	/*
	if (abs(i - n11) < abs(j - n11)) right = 1. - left + q;
	else left = 1.0 - right + q;
	*_left = left; *_right = right;
	return q;
	*/
}

void WorkMem::ReformatData(std::string &s, std::vector<int> &q, std::vector<int> &weight, std::vector<double> &condLike) {
	for (int i = 0; i < 4; ++i)
		QualityCount[i].clear();

	unsigned missingCount = 0;
	for(unsigned i = 0; i < s.length(); ++i) {
		int quality = q[i];
		if (s[i] == '?'){
			++missingCount;
		}else{
			++QualityCount[convert_table[s[i]]][quality];
		}
	}
	unsigned totalPattern = unsigned(QualityCount[0].size() + QualityCount[1].size() + QualityCount[2].size() + QualityCount[3].size());
	if(missingCount != 0)
		totalPattern += 1;

	for (int i = 0; i < 4; ++i){
		for (auto& x: QualityCount[i]){
			weight.push_back(x.second);
			double probIncorrect = pow(10.0, (double)(x.first)/(-10.0));
			auto push_val = probIncorrect/3.0;
			for (int j = 0; j < 4; ++j)
				condLike.push_back(push_val);

			condLike[condLike.size() - 4 + i] = 1.0 - probIncorrect;
		}
	}
	if (missingCount){
		weight.push_back(missingCount);
		for (int i = 0; i < 4; ++i)
			condLike.push_back(1.0);
	}

}

void WorkMem::ResetParams() {
	weight.clear();
	condLike.clear();
	pMatrix.assign(16, 0.0);
	baseFreq.assign(4, 0.25);
	kappa  = 1.0;
	brlens = BRLENS_PR_MEAN;
	lnlike = 0.0;

	unscaledBaseFreq.assign(4, 25.0);
	brlensMin = 1.E-8;
	brlensMax = 100.0;
	min_diff_lk_local = 1.E-04;
	optimizeBaseFreq  = false;
	
}

void WorkMem::Check_Br_Len_Bounds() {
	if(brlens > brlensMax) brlens = brlensMax;
	if(brlens < brlensMin) brlens = brlensMin;
}

void WorkMem::NormalizeBaseFreq() {
	double sum;
	int i;
	
	if(optimizeBaseFreq) {
		sum = .0;
		For(i,4) sum += fabs(unscaledBaseFreq[i]);
		For(i,4) baseFreq[i] = fabs(unscaledBaseFreq[i])/sum;
		
		do {
			sum = .0;
			For(i,4) {
				if(baseFreq[i] < 0.000001) baseFreq[i]=0.000001;
				if(baseFreq[i] > 0.999999) baseFreq[i]=0.999999;
				sum += baseFreq[i];
			}
			For(i,4) baseFreq[i]/=sum;
		}
		while((sum > 1.000001) || (sum < 0.999999));
	}
}

void WorkMem::CalcProbMat(std::vector<double> &pMat, const std::vector<double> &state_freqs, const double &kappa, const double &brlens) {
	pMat.assign(16, 0.0);
	
	double piA = state_freqs[0];
	double piC = state_freqs[1];
	double piG = state_freqs[2];
	double piT = state_freqs[3];
	
	double PiA = piA + piG;
	double PiC = piC + piT;
	double PiG = piA + piG;
	double PiT = piC + piT;
	
	double bigPiInvA = 1.0/PiA;
	double bigPiInvC = 1.0/PiC;
	double bigPiInvG = 1.0/PiG;
	double bigPiInvT = 1.0/PiT;
	
	double ta, tb, tc, td, y;
	double denom = ((piA + piG)*(piC + piT) + kappa*((piA*piG) + (piC*piT)));
	double beta  = 0.5/denom;
	
	// transition probability given the branch length, kappa, frequencies
	//
	double t = 0.5*brlens;
    // The next two lines fix the "Rota" bug; see BUGS file for details
//    if(t < 1.e-8)
//        t = 1.e-8; //TreeNode::edgeLenEpsilon;
	double x = exp(-beta*t);
	
	// changes to base A
	td			= -beta*(1 + PiA*(kappa - 1.0));
	y			= exp(t*td);
	ta			= piA*(bigPiInvA - 1.0);
	tb			= (PiA - piA)*bigPiInvA;
	tc			= piA*bigPiInvA;
	pMat[0*4+0]	= piA + (x*ta) + (y*tb);
	pMat[1*4+0]	= piA*(1.0 - x);
	pMat[2*4+0]	= piA + (x*ta) - (y*tc);
	pMat[3*4+0]	= pMat[1*4+0];
	
	// changes to base C
	td			= -beta*(1 + PiC*(kappa - 1.0));
	y			= exp(t*td);
	ta			= piC*(bigPiInvC - 1.0);
	tb			= (PiC - piC)*bigPiInvC;
	tc			= piC*bigPiInvC;
	pMat[0*4+1] = piC*(1.0 - x);
	pMat[1*4+1] = piC + (x*ta) + (y*tb);
	pMat[2*4+1] = pMat[0*4+1];
	pMat[3*4+1] = piC + (x*ta) - (y*tc);
	
	// changes to base G
	td			= -beta*(1 + PiG*(kappa - 1.0));
	y			= exp(t*td);
	ta			= piG*(bigPiInvG - 1.0);
	tb			= (PiG - piG)*bigPiInvG;
	tc			= piG*bigPiInvG;
	pMat[0*4+2] = piG + (x*ta) - (y*tc);
	pMat[1*4+2] = piG*(1.0 - x);
	pMat[2*4+2] = piG + (x*ta) + (y*tb);
	pMat[3*4+2] = pMat[1*4+2];
	
	// changes to base T
	td			= -beta*(1 + PiT*(kappa - 1.0));
	y			= exp(t*td);
	ta			= piT*(bigPiInvT - 1.0);
	tb			= (PiT - piT)*bigPiInvT;
	tc			= piT*bigPiInvT;
	pMat[0*4+3] = piT*(1.0 - x);
	pMat[1*4+3] = piT + (x*ta) - (y*tc);
	pMat[2*4+3] = pMat[0*4+3];
	pMat[3*4+3] = piT + (x*ta) + (y*tb);
}

void WorkMem::CalcLnLike(double &lnlike, const std::vector<int> &weight, const std::vector<double> &condLike, 
						const std::vector<double> &baseFreq, const std::vector<double> &pMatrix) {
	unsigned totalPattern = unsigned(weight.size());
	std::vector<double>::const_iterator prob  = pMatrix.begin();
	std::vector<double>::const_iterator input = condLike.begin();
	std::vector<int>::const_iterator weightIt = weight.begin();
	double tmpLnlike = 0.0;
	
	for(unsigned i = 0; i < totalPattern; i++) {
		double rootCondLike[4] = {0.0, 0.0, 0.0, 0.0};
		rootCondLike[0] = baseFreq[0] * ((*(prob   ))*(*(input)) + (*(prob+1 ))*(*(input+1)) + (*(prob+2 ))*(*(input+2)) + (*(prob+3 ))*(*(input+3)));
		rootCondLike[1] = baseFreq[1] * ((*(prob+4 ))*(*(input)) + (*(prob+5 ))*(*(input+1)) + (*(prob+6 ))*(*(input+2)) + (*(prob+7 ))*(*(input+3)));
		rootCondLike[2] = baseFreq[2] * ((*(prob+8 ))*(*(input)) + (*(prob+9 ))*(*(input+1)) + (*(prob+10))*(*(input+2)) + (*(prob+11))*(*(input+3)));
		rootCondLike[3] = baseFreq[3] * ((*(prob+12))*(*(input)) + (*(prob+13))*(*(input+1)) + (*(prob+14))*(*(input+2)) + (*(prob+15))*(*(input+3)));
		if(refLetter=='A') {
			rootCondLike[0] *= (*(prob   ));
			rootCondLike[1] *= (*(prob+4 ));
			rootCondLike[2] *= (*(prob+8 ));
			rootCondLike[3] *= (*(prob+12));
		}
		else if(refLetter=='C') {
			rootCondLike[0] *= (*(prob+1 ));
			rootCondLike[1] *= (*(prob+5 ));
			rootCondLike[2] *= (*(prob+9 ));
			rootCondLike[3] *= (*(prob+13));
		}
		else if(refLetter=='G') {
			rootCondLike[0] *= (*(prob+2 ));
			rootCondLike[1] *= (*(prob+6 ));
			rootCondLike[2] *= (*(prob+10));
			rootCondLike[3] *= (*(prob+14));
		}
		else if(refLetter=='T') {
			rootCondLike[0] *= (*(prob+3 ));
			rootCondLike[1] *= (*(prob+7 ));
			rootCondLike[2] *= (*(prob+11));
			rootCondLike[3] *= (*(prob+15));
		}
		tmpLnlike += (double)(*weightIt) * log(rootCondLike[0]+rootCondLike[1]+rootCondLike[2]+rootCondLike[3]);
		weightIt++;
		input += 4;
	}
	lnlike = tmpLnlike;

//#if POSTERIORMODE
//#if DREAM
	lnlike += (log(6.0) - brlens/BRLENS_PR_MEAN - log(BRLENS_PR_MEAN));
//#else
//	lnlike += (log(6.0) - kappa - brlens/BRLENS_PR_MEAN - log(BRLENS_PR_MEAN));
//#endif
//#endif
}

void WorkMem::Adjust_Min_Diff_Lk() {
	int exponent;
	exponent = (int)floor(log10(fabs(lnlike)));
	if(sizeof(double) == 8) {
		min_diff_lk_local = pow(10.,exponent - DBL_DIG + 1);
    }else if(sizeof(double) == 4) {
		min_diff_lk_local = pow(10.,exponent - FLT_DIG + 1);
    }
}

double WorkMem::Lk() {
	Check_Br_Len_Bounds();
	NormalizeBaseFreq();
	CalcProbMat(pMatrix, baseFreq, kappa, brlens);
	CalcLnLike(lnlike, weight, condLike, baseFreq, pMatrix);
	Adjust_Min_Diff_Lk();
	return lnlike;
}

double WorkMem::Generic_Brent_Lk(double *param, double ax, double cx, double tol, int n_iter_max, int quickdirty) {
	int    iter;
	double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
	double e = 0.0;
	double old_lnL, init_lnL;
	double bx = *param;
	
	d = 0.0;
	a = ((ax < cx) ? ax : cx);
	b = ((ax > cx) ? ax : cx);
	x = w = v = bx;
	old_lnL  = UNLIKELY;
	(*param) = bx;
	fw = fv = fx = fu = -Lk();
	init_lnL = -fw;
	
	for(iter=1; iter<=BRENT_IT_MAX; iter++) {
		xm   = 0.5*(a+b);
		tol2 = 2.0*(tol1=tol*x+BRENT_ZEPS);
		
		if((fu > init_lnL + tol) && (quickdirty)) {
			(*param) = x;
			fu = Lk();
			return fu;
		}
		
		if((fabs(fu-old_lnL) < tol) || (iter > n_iter_max - 1)) {
			(*param) = x;
			fu = Lk();
			return fu;
		}
		
		if(fabs(e) > tol1) {
			r = (x-w)*(fx-fv);
			q = (x-v)*(fx-fw);
			p = (x-v)*q-(x-w)*r;
			q = 2.0*(q-r);
			if(q > 0.0)
				p = -p;
			q     = fabs(q);
			etemp = e;
			e     = d;
			if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)) {
				d = BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
			}
			else {
				d = p/q;
				u = x+d;
				if(u-a < tol2 || b-u < tol2)
					d = SIGN(tol1,xm-x);
				/* PhyML_Printf(" Parabolic step [e=%f]\n",e); */
			}
        }
		else {
			d = BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
			/* PhyML_Printf(" Golden section step (default) [e=%f tol1=%f a=%f b=%f d=%f]\n",e,tol1,a,b,d); */
		}
		
		u = (fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		(*param) = u;
		old_lnL  = fu;
		fu       = -Lk();
		
		if(fu <= fx) {
			if(u >= x)
				a = x;
			else
				b = x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		}
		else {
			if(u < x)
				a = u;
			else
				b = u;
			if(fu < fw || fabs(w-x) < SMALL) {
				v  = w;
				w  = u;
				fv = fw;
				fw = fu;
			}
			else if(fu < fv || fabs(v-x) < SMALL || fabs(v-w) < SMALL) {
				v  = u;
				fv = fu;
			}
		}
    }
	
	std::cerr << "\nToo many iterations in BRENT!\n";
	exit(1);
	return(-1);
}

double WorkMem::Br_Len_Brent(double *brlens) {
	Generic_Brent_Lk(brlens,
					 brlensMin,
					 brlensMax,
					 min_diff_lk_local,
					 BRENT_IT_MAX,
					 0);
	return lnlike;
}

void WorkMem::Optimize_Br_Len_Serie() {
	double lk_init = lnlike;
	Br_Len_Brent(&brlens);
	if(lnlike < lk_init - min_diff_lk_local) {
		std::cerr << "\nError in Optimize_Br_Len_Serie.\n";
		exit(1);
    }
}

void WorkMem::Optimize_State_Freqs() {
	optimizeBaseFreq = true;

	int i;
	For(i,4) {
		Generic_Brent_Lk(&(unscaledBaseFreq[i]),
						 0.,
						 100.,
						 min_diff_lk_local,
						 BRENT_IT_MAX,
						 0);
	}
	
	optimizeBaseFreq = false;
}

void WorkMem::Optimiz_All_Free_Param() {
//#if !DREAM
//	Optimize_TsTv();
//#endif
	Optimize_State_Freqs();
	Lk();
}


void WorkMem::Round_Optimize(int n_round_max) {
	int    n_round, each;
	double lk_old, lk_new;
	
	lk_new  = UNLIKELY;
	lk_old  = UNLIKELY;
	n_round = 0;
	each    = 0;
	
	Lk();
	
	while(n_round < n_round_max) {
		Optimize_Br_Len_Serie();
		Lk();
		
		if(!each) {
			each = 1;
			Optimiz_All_Free_Param();
		}
		
		lk_new = lnlike;
//		if(lk_new < lk_old - min_diff_lk_local) {
//			std::cerr << lk_new << " " << lk_old << " " << min_diff_lk_local << "\nOptimisation failed ! (Round_Optimize).\n";
//			exit(1);
//		}
		if(fabs(lk_new - lk_old) < min_diff_lk_local)
			break;
		else
			lk_old = lk_new;
		n_round++;
		each--;
    }
	
	Optimiz_All_Free_Param();
}

void WorkMem::CalcEvolDistance(std::string &s, std::vector<int> &q) {
	ResetParams();
	ReformatData(s, q, weight, condLike);
	Round_Optimize(ROUND_MAX);
}

string WorkMem::scientificStr(double num){
	os.str("");
	os.clear();
	os << num;
	return os.str();
}


#define bam1_aux(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname + (b)->core.l_qseq + ((b)->core.l_qseq + 1)/2)

inline int bam_aux_type2size(int x) {
	if (x == 'C' || x == 'c' || x == 'A') return 1;
	else if (x == 'S' || x == 's') return 2;
	else if (x == 'I' || x == 'i' || x == 'f') return 4;
	else return 0;
}

#define __skip_tag(s) do { \
int type = (*(s)); \
++(s); \
if (type == 'Z' || type == 'H') { while (*(s)) ++(s); ++(s); } \
else if (type == 'B') (s) += 5 + bam_aux_type2size(*(s)) * (*(int32_t*)((s)+1)); \
else (s) += bam_aux_type2size(type); \
} while(0)

uint8_t *bam_aux_get_muse(const bam1_t *b, const char tag[2]) {
	uint8_t *s;
	int y = tag[0]<<8 | tag[1];
	s = bam1_aux(b);
	while(s < b->data + b->l_data) {
		int x = (int)s[0]<<8 | s[1];
		s += 2;
		if(x == y)
			return s;
		__skip_tag(s);
	}
	return 0;
}

void WorkMem::processPileupFun(bam_hdr_t* h, mplp_conf_t* conf, std::string& output, int tid, uint32_t pos, const char* ref_ptr, uint32_t ref_len){

	refLetter = ref_ptr[pos];
	if (refLetter == 'N' || refLetter == '?')
		return;
	altLetter = ' ';

	for (int i = 0; i < plp.size(); ++i){
		for (int j = 0; j < plp[i].size(); ++j)
			resolve_cigar2_pb2(plp[i][j], pos, plp[i][j]->s);
	}

	rawS[0].clear();rawS[1].clear();
	rawQ[0].clear();rawQ[1].clear();
	overlapPair[0].clear(); overlapPair[1].clear();
	bool ifAlt = false;
	for (int i = 0; i < 2; ++i){
		for(int j = 0; j < plp[i].size(); ++j) {
			auto p = plp[i][j];
			bool ifAlt_tmp = pileup_seq(p, refLetter, rawS[i]);
			if (i == 0)
				ifAlt |= ifAlt_tmp;
			rawQ[i].push_back(min(93, int(bam_get_qual(p->b)[p->qpos])));
			std::pair<const char*, size_t> key = make_pair((char*)(p->b->data), charArrayHashFun((char*)(p->b->data)));

			auto iter = overlapPair[i].find(key);
			if (iter != overlapPair[i].end()){
				iter->second.second = j;
			}else{
				overlapPair[i][move(key)] = make_pair(j, -1);
			}
		}
	}

	if (!ifAlt)
		return;

	rawTumorPileupElementPointer = plp[0];
	
	for (int i = 0; i < 2; ++i){
		for (auto& one_pair: overlapPair[i]){
			if (one_pair.second.second != -1){
				if(rawS[i][(one_pair.second).first] != rawS[i][(one_pair.second).second]) {
					if (i == 0){
						plp[i][one_pair.second.first] = nullptr;
						plp[i][one_pair.second.second] = nullptr;
					}else{
						if(rawS[i][(one_pair.second).first] == refLetter){
							plp[i][(one_pair.second).first] = nullptr;
						}else if(rawS[i][(one_pair.second).second] == refLetter){
							plp[i][(one_pair.second).second] = nullptr;
						}
					}
				}else{
					if(rawQ[i][(one_pair.second).first] < rawQ[i][(one_pair.second).second]){
						plp[i][(one_pair.second).first] = nullptr;
					}else{
						plp[i][(one_pair.second).second] = nullptr;
					}
				}
			}
		}
		//ignore here
	}

	// proximal gap filter
	//
	int  deletionCount  = 0;
	int  insertionCount = 0;

	for (auto one_plp: plp[0]){
		if (!one_plp) continue;
		if(deletionCount >= GAP_EVENT_CUTOFF || insertionCount >= GAP_EVENT_CUTOFF)
			break;
		if(one_plp->is_del) {
			++deletionCount;
		}else{
			uint32_t *cigar = bam_get_cigar(one_plp->b);
			int positionInRead = one_plp->qpos;
			int operationStart = 0;
			for(int j = 0; j < one_plp->b->core.n_cigar; ++j) {
				int operation = _cop(cigar[j]);
				int length    = _cln(cigar[j]);
				if(operation==BAM_CINS) {
					int distance = (operationStart>positionInRead) ? operationStart-positionInRead : positionInRead-(operationStart+length-1);
					if(distance<=GAP_EVENT_PROXIMITY) {
						++insertionCount;
						break;
					}
				}

				if(operation==BAM_CDEL) {
					int distance = (operationStart>positionInRead) ? operationStart-positionInRead : positionInRead-operationStart+1;
					if(distance<=GAP_EVENT_PROXIMITY) {
						++deletionCount;
						break;
					}
				}else
					operationStart += length;
			}
		}
	}

	if(deletionCount>=GAP_EVENT_CUTOFF || insertionCount>=GAP_EVENT_CUTOFF)
		return;

	s[0].clear(); s[1].clear();
	q[0].clear(); q[1].clear();
	noOverlapNoDelPassMinQNoMapQ0TumorPileupElementPointer.clear();
	for(int i = 0; i < plp[0].size(); ++i) {
		if (!plp[0][i])
			continue;
		auto p = plp[0][i];
		if((rawS[0][i]!='?') && (!(p->is_del)) && (p->b->core.flag&1) && (!(p->b->core.flag&8)) && (p->b->core.qual > 0) && (rawQ[0][i] >= MIN_QUALITY_SCORE)) {
			noOverlapNoDelPassMinQNoMapQ0TumorPileupElementPointer.push_back(p);
			s[0].push_back(rawS[0][i]);
			q[0].push_back(rawQ[0][i]);
		}
	}
	finalNormalPileupElementPointer.clear();
	for(int i = 0; i < plp[1].size(); ++i) {
		if (!plp[1][i])
			continue;
		auto p = plp[1][i];
		if((rawS[1][i]!='?') && (!(p->is_del)) && (rawQ[1][i] >= MIN_QUALITY_SCORE)) {
			finalNormalPileupElementPointer.push_back(p);
			s[1] += rawS[1][i];
			q[1].push_back(rawQ[1][i]);
		}
	}

	if(s[0].length() < conf->depth_cutoff || s[1].length() < conf->depth_cutoff) 
		return;

	for (int i = 0; i < 4; ++i){
		tumorVariantAlleleCount[i] = std::make_pair(i, 0);
	}
	double obsTVAF = -1.0;
	double obsNVAF = -1.0;
	int nRefInTumor = 0;
	ifAlt = false;
	for(int i = 0; i < s[0].length(); i++) {
		if(s[0][i]==refLetter) {
			nRefInTumor += 1;
		}
		else if(s[0][i]!='?') {
			tumorVariantAlleleCount[convert_table[s[0][i]]].second += 1;
			ifAlt = true;
		}
	}

	if (!ifAlt) return;

	sort(tumorVariantAlleleCount, tumorVariantAlleleCount + 4, 
		[](const std::pair<int, int>& left, const std::pair<int, int>& right)->bool{return left.second > right.second;});
	
	bool useSecondMutantAllele = false;
	int nRefInNormal = 0;

	for (int i = 0; i < 4; ++i)
		normalVariantAlleleCountAndQuality[i] = std::make_pair(0, 0);

	for(int i = 0; i < s[1].length(); ++i) {
		if(s[1][i]==refLetter) {
			nRefInNormal += 1;
		}else if(s[1][i]!='?') {
			auto& tmp_pair = normalVariantAlleleCountAndQuality[convert_table[s[1][i]]];
			tmp_pair.first   += 1;
			tmp_pair.second  += q[1][i];
		}
	}

	auto altLetterIdx = tumorVariantAlleleCount[0].first;
	altLetter = revert_table[altLetterIdx];
	obsTVAF = double(tumorVariantAlleleCount[0].second)/double(s[0].length());
	obsNVAF = double(normalVariantAlleleCountAndQuality[altLetterIdx].first)/double(s[1].length());

	bool observedInNormal = false;
	
	if((normalVariantAlleleCountAndQuality[altLetterIdx].first >=NORMAL_OBSERVATION_COUNT || obsNVAF>=NORMAL_OBSERVATION_VAF) &&
		normalVariantAlleleCountAndQuality[altLetterIdx].second > NORMAL_OBSERVATION_SUM_QUAL){
		observedInNormal = true;
	}

	if(observedInNormal && tumorVariantAlleleCount[1].second != 0) {
		auto secondAlleleIdx = tumorVariantAlleleCount[1].first;
		obsTVAF = double(tumorVariantAlleleCount[1].second)/double(s[0].length());
		obsNVAF = double(normalVariantAlleleCountAndQuality[secondAlleleIdx].first)/double(s[1].length());
		if(obsTVAF >= 0.05) {
			if(!((normalVariantAlleleCountAndQuality[secondAlleleIdx].first>=NORMAL_OBSERVATION_COUNT || obsNVAF>=NORMAL_OBSERVATION_VAF) && 
				  normalVariantAlleleCountAndQuality[secondAlleleIdx].second>NORMAL_OBSERVATION_SUM_QUAL)) {

				double hetTwotail, homTwotail;
				int coverageWithoutSecondAltAllele = s[1].length() - normalVariantAlleleCountAndQuality[secondAlleleIdx].first;


				kt_fisher_exact(coverageWithoutSecondAltAllele/2, coverageWithoutSecondAltAllele/2, 
								normalVariantAlleleCountAndQuality[altLetterIdx].first, 
								coverageWithoutSecondAltAllele - normalVariantAlleleCountAndQuality[altLetterIdx].first,
								&hetTwotail);


				kt_fisher_exact(coverageWithoutSecondAltAllele, 0, normalVariantAlleleCountAndQuality[altLetterIdx].first, 
								coverageWithoutSecondAltAllele - normalVariantAlleleCountAndQuality[altLetterIdx].first, &homTwotail);

				if(hetTwotail >= GERMLINE_VARIANT_FET_CUTOFF || homTwotail >= GERMLINE_VARIANT_FET_CUTOFF) {
					observedInNormal = false;
					useSecondMutantAllele = true;
					altLetter = revert_table[secondAlleleIdx];
				}

			}
		}
	}
	
	if (observedInNormal)
		return;
	
	if(obsTVAF < conf->minAltFraction)
		return;

	int strandBiasContingencyTable[2][2] = {{0,0}, {0,0}};
	for(int i = 0; i < rawS[0].length(); ++i) {
		auto p = rawTumorPileupElementPointer[i];
		if((rawS[0][i]!='?') && (!(p->is_del)) && (p->b->core.flag&1) && (!(p->b->core.flag&8)) && (p->b->core.qual > 0) && (rawQ[0][i] >= MIN_QUALITY_SCORE)) {
			if(rawS[0][i] != altLetter) {
				if(bam_is_rev(p->b))
					strandBiasContingencyTable[1][1] += 1;
				else
					strandBiasContingencyTable[1][0] += 1;
			}
			else {
				if(bam_is_rev(p->b))
					strandBiasContingencyTable[0][1] += 1;
				else
					strandBiasContingencyTable[0][0] += 1;
			}
		}
	}
	double twotail;
	kt_fisher_exact(strandBiasContingencyTable[0][0], strandBiasContingencyTable[0][1], strandBiasContingencyTable[1][0], strandBiasContingencyTable[1][1], &twotail);
	if(twotail <= STRAND_BIAS_FET_CUTOFF)
		return;
	
	auto& finalTumorPileupElementPointer = rawTumorPileupElementPointer;
	finalTumorPileupElementPointer.clear();

	std::string tmpTumorS;
	tmpTumorS.swap(s[0]);
	//s[0] = "";
	std::vector<int> tmpTumorQ;
	tmpTumorQ.swap(q[0]);
	//q[0].clear();
	//q[0].reserve(n_plp[0]);

	bool existRequiredRead = false;
	std::pair<int, int>& tumorSomaticMutantReadMAPQ = tumorVariantAlleleCount[0];
	tumorSomaticMutantReadMAPQ.first = tumorSomaticMutantReadMAPQ.second = 0;

	for (int i = 0; i < noOverlapNoDelPassMinQNoMapQ0TumorPileupElementPointer.size(); ++i){
		bam_pileup1_t_pb *p = noOverlapNoDelPassMinQNoMapQ0TumorPileupElementPointer[i];
		if(tmpTumorS[i] == altLetter) {
			//tumorSomaticMutantReadMAPQ.push_back(mapQ);
			tumorSomaticMutantReadMAPQ.first += p->b->core.qual;
			++tumorSomaticMutantReadMAPQ.second;
		}
	}

	if ((double(tumorSomaticMutantReadMAPQ.first) / tumorSomaticMutantReadMAPQ.second) <= LOW_MEAN_MAPQ)
		return;

	smallestDistance.clear();
	
	for(int i = 0; i < noOverlapNoDelPassMinQNoMapQ0TumorPileupElementPointer.size(); ++i) {
		bam_pileup1_t_pb *p = noOverlapNoDelPassMinQNoMapQ0TumorPileupElementPointer[i];

		if(p->b->core.qual<=LOW_MAPQ)
			continue;
		uint8_t *symbol = bam_aux_get_muse(p->b, "XT");
		if(symbol) {
			if(char(*(symbol+1))=='M')
				continue;
		}

		uint32_t softclipLength = 0;
		uint32_t *cigar  = bam_get_cigar(p->b);
		for(int ithCigar = 0; ithCigar < p->b->core.n_cigar; ++ithCigar) {
			int op = cigar[ithCigar] & BAM_CIGAR_MASK;
			if(op == BAM_CSOFT_CLIP)
				softclipLength += cigar[ithCigar] >> BAM_CIGAR_SHIFT;
		}
		if((double)softclipLength/(double)(p->b->core.l_qseq) >= HEAVY_SOFTCLIP)
			continue;
		uint8_t     *seq   = bam_get_seq(p->b);
		uint8_t     *qual  = bam_get_qual(p->b);
		auto& core  = p->b->core;
		int mismatchSumQual = 0;
		for(int k = 0, y = 0, x = core.pos; k < core.n_cigar; ++k) {
			int j, l = _cln(cigar[k]), op = _cop(cigar[k]);
			if(op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
				for(j = 0; j < l; ++j) {
					int z = y + j;
					int c1 = bam_seqi(seq, z); 
					//if((*ref_ptr)[x+j] == 0)
					if (x + j >= ref_len)
						break; // out of boundary
					int c2 = bam_nt16_table[(ref_ptr)[x+j]];
					if((c1 == c2 && c1 != 15 && c2 != 15) || c1 == 0) { // a match
						continue;
					}else {
						mismatchSumQual += int(qual[z]);
					}
				}
				if(j < l)
					break;
				x += l;
				y += l;
			}
			else if(op == BAM_CDEL) {
				/*
				for(j = 0; j < l; ++j) {
					if((*ref_ptr)[x+j] == 0)
						break;
				}
				if(j < l)
					break;
				*/
				if (x + l - 1 >= ref_len)
					break;
				x += l;
			}
			else if(op == BAM_CINS || op == BAM_CSOFT_CLIP) {
				y += l;
			}
			else if(op == BAM_CREF_SKIP) {
				x += l;
			}
		}
		if(mismatchSumQual > HEAVY_MISMATCH)
			continue;

		finalTumorPileupElementPointer.push_back(p);
		s[0].push_back(tmpTumorS[i]);
		q[0].push_back(tmpTumorQ[i]);

		if(tmpTumorS[i]==altLetter){
			int minDistance = min(abs(int(pos) - int(core.pos)), int(abs(int64_t(pos) - p->end)));
			smallestDistance.push_back(minDistance);
			if(symbol) {
				if(char(*(symbol+1))=='R')
					continue;
			}

			if((core.flag&1) && (core.flag&2)) {
				if(core.qual >= REQUIRED_MAPQ && tmpTumorQ[i] >= REQUIRED_BASE_QUALITY) {
					existRequiredRead = true;
				}
			}
		}
	}

	if (!existRequiredRead)
		return;

	auto nMutant = smallestDistance.size();
	std::partial_sort(smallestDistance.begin(), smallestDistance.begin() + (nMutant / 2) + 1, smallestDistance.end());
	double distanceMedian = (nMutant % 2 == 1) ? smallestDistance[nMutant / 2]:
												(smallestDistance[(nMutant / 2) - 1] + smallestDistance[nMutant / 2]) / 2;
	for (auto& x: smallestDistance){
		x = fabs(x - distanceMedian);
	}
	std::partial_sort(smallestDistance.begin(), smallestDistance.begin() + (nMutant / 2) + 1, smallestDistance.end());
	double distanceMAD = (nMutant % 2 == 1) ? smallestDistance[nMutant / 2]:
											(smallestDistance[(nMutant / 2) - 1] + smallestDistance[nMutant / 2]) / 2;

	if(distanceMedian <= CLUSTERED_POSITION_MEDIAN_CUTOFF && distanceMAD<=CLUSTERED_POSITION_MAD_CUTOFF)
		return;

	overlapRG.clear();
	for(int i = 0; i < s[0].length(); i++) {
		if(s[0][i]==altLetter) {
			const uint8_t *rg = bam_aux_get_muse(finalTumorPileupElementPointer[i]->b, "RG");
			if(rg) {
				overlapRG.insert(string((char*)(rg+1)));
			}
		}
	}
	int nGoodRefInNormal = 0;
	for(int i = 0; i < s[1].length(); ++i) {
		if(s[1][i]==refLetter) {
			// Mapping Quality score > 0
			//
			auto& core = finalNormalPileupElementPointer[i]->b->core;
			if (core.qual > 0 && ((core.flag & 3) == 3) && (!(core.flag & 8)))
				++nGoodRefInNormal;
		}
	}

	CalcEvolDistance(s[0], q[0]);

	if(brlens <= conf->min_output_brlens)
		return;

	if(std::count(s[1].begin(), s[1].end(), altLetter) != 0) {
		std::vector<double>& tumorBaseFreq = smallestDistance;
		tumorBaseFreq.swap(baseFreq);

		double tumorKappa  = kappa;
		double tumorBrlens = brlens;
		CalcEvolDistance(s[1], q[1]);
		output += string(h->target_name[tid]) + "\t" + to_string(pos + 1) + "\t" + refLetter + "\t" + altLetter + "\t" + to_string(s[0].length()) +
					 "\t" + to_string(s[1].length()) + "\t";
		output += /*std::scientific +*/ scientificStr(obsTVAF) + "\t" + scientificStr(tumorBaseFreq[0]) + "\t" + scientificStr(tumorBaseFreq[1]) + "\t" + scientificStr(tumorBaseFreq[2]) + 
					 			"\t" + scientificStr(tumorBaseFreq[3]) + "\t" + scientificStr(tumorKappa) + "\t" + scientificStr(tumorBrlens) + "\t";
		output += /*std::scientific +*/ scientificStr(obsNVAF) + "\t" + scientificStr(baseFreq[0]) + "\t" + scientificStr(baseFreq[1]) + "\t" + scientificStr(baseFreq[2]) + 
					 			"\t" + scientificStr(baseFreq[3]) + "\t" + scientificStr(kappa) + "\t" + scientificStr(brlens) + "\t" + to_string(nGoodRefInNormal) + "\t" + to_string(overlapRG.size()) + "\t";

		if(useSecondMutantAllele)
			output += string("Y\t");
		else
			output += string("N\t");
		output += to_string(std::count(s[0].begin(), s[0].end(), refLetter)) + "\t" + to_string(std::count(s[0].begin(), s[0].end(), altLetter)) + "\t" +
					to_string(std::count(s[1].begin(), s[1].end(), refLetter)) + "\t" + to_string(std::count(s[1].begin(), s[1].end(), altLetter)) + "\t";
	
	}else{
		output += string(h->target_name[tid]) + "\t" + to_string(pos + 1) + "\t" + refLetter + "\t" + altLetter + "\t" + to_string(s[0].length()) + "\t" + to_string(s[1].length()) + "\t";
		output += /*std::scientific <<*/ scientificStr(obsTVAF) + "\t" + scientificStr(baseFreq[0]) + "\t" + scientificStr(baseFreq[1]) + "\t" + scientificStr(baseFreq[2]) + "\t" + scientificStr(baseFreq[3]) + "\t"
				 + scientificStr(kappa) + "\t" + scientificStr(brlens) + "\t";
		output += /*std::scientific <<*/ scientificStr(obsNVAF) + "\t";
		output += string("NA\tNA\tNA\tNA\tNA\t") + /*to_string(MIN_BRLENS)*/ "1.00000e-08" + "\t" + to_string(nGoodRefInNormal) + "\t" + to_string(overlapRG.size()) + "\t";
		if(useSecondMutantAllele)
			output += "Y\t";
		else
			output += "N\t";
		output += to_string(std::count(s[0].begin(), s[0].end(), refLetter)) + "\t" + to_string(std::count(s[0].begin(), s[0].end(), altLetter)) + "\t" +
				to_string(std::count(s[1].begin(), s[1].end(), refLetter)) + "\t" + to_string(std::count(s[1].begin(), s[1].end(), altLetter)) + "\t";

	}

	altString.clear();
	normalFormat.clear();
	tumorFormat.clear();

	altString += altLetter;
	for (int i = 0; i < 5; ++i){
		rawNormalAlleleCount[i] = 0;
		rawNormalAlleleQuality[i] = 0;
		rawTumorAlleleCount[i] = 0;
		rawTumorAlleleQuality[i] = 0;
	}

	for(int i = 0; i < rawS[1].length(); ++i) {
		rawNormalAlleleCount[convert_table[rawS[1][i]]]   += 1;
		rawNormalAlleleQuality[convert_table[rawS[1][i]]] += rawQ[1][i];
	}

	for(int i = 0; i < rawS[0].length(); ++i) {
		rawTumorAlleleCount[convert_table[rawS[0][i]]]   += 1;
		rawTumorAlleleQuality[convert_table[rawS[0][i]]] += rawQ[0][i];
	}

	if((rawNormalAlleleCount[convert_table[refLetter]]+rawNormalAlleleCount[convert_table[altLetter]]+rawNormalAlleleCount[convert_table['?']]) == rawS[1].length()) {
		normalFormat = "0/0:" + to_string(rawS[1].length()-rawNormalAlleleCount[convert_table['?']]) + ":" + to_string(rawNormalAlleleCount[convert_table[refLetter]]) + "," + 
								to_string(rawNormalAlleleCount[convert_table[altLetter]]) + ":";
		if(rawNormalAlleleCount[convert_table[refLetter]]==0) {
			normalFormat += "0,";
		}
		else {
			int tmp = round((double)rawNormalAlleleQuality[convert_table[refLetter]]/(double)rawNormalAlleleCount[convert_table[refLetter]]);
			normalFormat += (to_string(tmp) + ",");
		}
		if(rawNormalAlleleCount[convert_table[altLetter]]==0) {
			normalFormat += "0:.";
		}
		else {
			int tmp = round((double)rawNormalAlleleQuality[convert_table[altLetter]]/(double)rawNormalAlleleCount[convert_table[altLetter]]);
			normalFormat += (to_string(tmp) + ":.");
		}

		int coverage = rawTumorAlleleCount[convert_table[refLetter]] + rawTumorAlleleCount[convert_table[altLetter]];

		double homTwotail;
		kt_fisher_exact(coverage, 0, rawTumorAlleleCount[convert_table[altLetter]], rawTumorAlleleCount[convert_table[refLetter]], &homTwotail);
		if(homTwotail>=GERMLINE_VARIANT_FET_CUTOFF) {
			// the test result is not significant, so I cannot reject the hypothesis that this position is homo-variant
			//
			tumorFormat = "1/1:";
		}
		else {
			tumorFormat = "0/1:";
		}

		tumorFormat += (to_string(rawS[0].length()-rawTumorAlleleCount[convert_table['?']]) + ":" + to_string(rawTumorAlleleCount[convert_table[refLetter]]) + "," + 
						to_string(rawTumorAlleleCount[convert_table[altLetter]]) + ":");

		if(rawTumorAlleleCount[convert_table[refLetter]]==0) {
			tumorFormat += "0,";
		}
		else {
			int tmp = round((double)rawTumorAlleleQuality[convert_table[refLetter]]/(double)rawTumorAlleleCount[convert_table[refLetter]]);
			tumorFormat += (to_string(tmp) + ",");
		}
		if(rawTumorAlleleCount[convert_table[altLetter]]==0) {
			tumorFormat += "0:2";
		}
		else {
			int tmp = round((double)rawTumorAlleleQuality[convert_table[altLetter]]/(double)rawTumorAlleleCount[convert_table[altLetter]]);
			tumorFormat += (to_string(tmp) + ":2");
		}

	}else{
		
		for (int i = 0; i < 5; ++i)
			rawNormalAlleleCountVec[i] = make_pair(revert_table[i], rawNormalAlleleCount[i]);
		
		std::sort(rawNormalAlleleCountVec, rawNormalAlleleCountVec + 4, 
				[](const std::pair<int, int>& left, const std::pair<int, int>& right)->bool{return left.second > right.second;});
		char potentialGermlineAllele = ' ';
		for(int i = 0; i < 5; i++) {
			char tmp_char = rawNormalAlleleCountVec[i].first;
			if( tmp_char!=refLetter && tmp_char!=altLetter && tmp_char!='?' && rawNormalAlleleCountVec[i].second>rawNormalAlleleCount[convert_table[altLetter]]) {
				potentialGermlineAllele = tmp_char;
				break;
			}
		}
		if(potentialGermlineAllele==' ') {
			normalFormat = "0/0:" + to_string(rawS[1].length()-rawNormalAlleleCount[convert_table['?']]) + ":" + to_string(rawNormalAlleleCount[convert_table[refLetter]]) + "," + 
									to_string(rawNormalAlleleCount[convert_table[altLetter]]) + ":";
			if(rawNormalAlleleCount[convert_table[refLetter]]==0) {
				normalFormat += "0,";
			}
			else {
				int tmp = round((double)rawNormalAlleleQuality[convert_table[refLetter]]/(double)rawNormalAlleleCount[convert_table[refLetter]]);
				normalFormat += (to_string(tmp) + ",");
			}
			if(rawNormalAlleleCount[convert_table[altLetter]]==0) {
				normalFormat += "0:.";
			}
			else {
				int tmp = round((double)rawNormalAlleleQuality[convert_table[altLetter]]/(double)rawNormalAlleleCount[convert_table[altLetter]]);
				normalFormat += (to_string(tmp) + ":.");
			}
			double homTwotail;
			int coverage = rawTumorAlleleCount[convert_table[refLetter]] + rawTumorAlleleCount[convert_table[altLetter]];

			kt_fisher_exact(coverage, 0, rawTumorAlleleCount[convert_table[altLetter]], rawTumorAlleleCount[convert_table[refLetter]], &homTwotail);
			
			if(homTwotail>=GERMLINE_VARIANT_FET_CUTOFF) {
				// the test result is not significant, so I cannot reject the hypothesis that this position is homo-variant
				//
				tumorFormat = "1/1:";
			}
			else {
				tumorFormat = "0/1:";
			}
			tumorFormat += (to_string(rawS[0].length()-rawTumorAlleleCount[convert_table['?']]) + ":" + to_string(rawTumorAlleleCount[convert_table[refLetter]]) + "," + 
							to_string(rawTumorAlleleCount[convert_table[altLetter]]) + ":");
			
			if(rawTumorAlleleCount[convert_table[refLetter]]==0) {
				tumorFormat += "0,";
			}
			else {
				int tmp = round((double)rawTumorAlleleQuality[convert_table[refLetter]]/(double)rawTumorAlleleCount[convert_table[refLetter]]);
				tumorFormat += (to_string(tmp) + ",");
			}
			if(rawTumorAlleleCount[convert_table[altLetter]]==0) {
				tumorFormat += "0:2";
			}
			else {
				int tmp = round((double)rawTumorAlleleQuality[convert_table[altLetter]]/(double)rawTumorAlleleCount[convert_table[altLetter]]);
				tumorFormat += (to_string(tmp) + ":2");
			}
		}else{
			double hetTwotail;
			int coverage = rawNormalAlleleCount[convert_table[refLetter]] + rawNormalAlleleCount[convert_table[potentialGermlineAllele]];

			kt_fisher_exact(coverage/2, coverage/2, rawNormalAlleleCount[convert_table[potentialGermlineAllele]], rawNormalAlleleCount[convert_table[refLetter]], &hetTwotail);
			
			if(hetTwotail>=GERMLINE_VARIANT_FET_CUTOFF) {
				std::string tmpS_tmp(1, potentialGermlineAllele);
				altString += ("," + tmpS_tmp);
				
				// normal Format first
				//
				normalFormat = "0/2:" + to_string(rawS[1].length()-rawNormalAlleleCount[convert_table['?']]) + ":" + to_string(rawNormalAlleleCount[convert_table[refLetter]]) + "," + 
								to_string(rawNormalAlleleCount[convert_table[altLetter]]) + "," + to_string(rawNormalAlleleCount[convert_table[potentialGermlineAllele]]) + ":";
				if(rawNormalAlleleCount[convert_table[refLetter]]==0) {
					normalFormat += "0,";
				}
				else {
					int tmp = round((double)rawNormalAlleleQuality[convert_table[refLetter]]/(double)rawNormalAlleleCount[convert_table[refLetter]]);
					normalFormat += (to_string(tmp) + ",");
				}
				if(rawNormalAlleleCount[convert_table[altLetter]]==0) {
					normalFormat += "0,";
				}
				else {
					int tmp = round((double)rawNormalAlleleQuality[convert_table[altLetter]]/(double)rawNormalAlleleCount[convert_table[altLetter]]);
					normalFormat += (to_string(tmp) + ",");
				}
				if(rawNormalAlleleCount[convert_table[potentialGermlineAllele]]==0) {
					normalFormat += "0:.";
				}
				else {
					int tmp = round((double)rawNormalAlleleQuality[convert_table[potentialGermlineAllele]]/(double)rawNormalAlleleCount[convert_table[potentialGermlineAllele]]);
					normalFormat += (to_string(tmp) + ":.");
				}
				// tumor Format
				//
				// always assume the GT is 0/1 unless test result shows otherwise
				//
				// n11       n12         | null (homo-variant)
				// n21       n22         | observation
				//-----------------------+-------------------------------
				// variant   ref         | Total
				//
				double homTwotail;
				int coverage = rawTumorAlleleCount[convert_table[refLetter]] + rawTumorAlleleCount[convert_table[altLetter]];

				kt_fisher_exact(coverage, 0, rawTumorAlleleCount[convert_table[altLetter]], rawTumorAlleleCount[convert_table[refLetter]], &homTwotail);
				
				if(homTwotail>=GERMLINE_VARIANT_FET_CUTOFF) {
					// the test result is not significant, so I cannot reject the hypothesis that this position is homo-variant
					//
					if(rawTumorAlleleCount[convert_table[potentialGermlineAllele]]!=0) {
						tumorFormat = "1/1/2:";
					}
					else {
						tumorFormat = "1/1:";
					}
				}
				else {
					if(rawTumorAlleleCount[convert_table[potentialGermlineAllele]]!=0) {
						tumorFormat = "0/1/2:";
					}
					else {
						tumorFormat = "0/1:";
					}
				}
				tumorFormat += (to_string(rawS[0].length()-rawTumorAlleleCount[convert_table['?']]) + ":" + to_string(rawTumorAlleleCount[convert_table[refLetter]]) + "," + 
								to_string(rawTumorAlleleCount[convert_table[altLetter]]) + "," + to_string(rawTumorAlleleCount[convert_table[potentialGermlineAllele]]) + ":");
				
				if(rawTumorAlleleCount[convert_table[refLetter]]==0) {
					tumorFormat += "0,";
				}
				else {
					int tmp = round((double)rawTumorAlleleQuality[convert_table[refLetter]]/(double)rawTumorAlleleCount[convert_table[refLetter]]);
					tumorFormat += (to_string(tmp) + ",");
				}
				if(rawTumorAlleleCount[convert_table[altLetter]]==0) {
					tumorFormat += "0,";
				}
				else {
					int tmp = round((double)rawTumorAlleleQuality[convert_table[altLetter]]/(double)rawTumorAlleleCount[convert_table[altLetter]]);
					tumorFormat += (to_string(tmp) + ",");
				}
				if(rawTumorAlleleCount[convert_table[potentialGermlineAllele]]==0) {
					tumorFormat += "0:2";
				}
				else {
					int tmp = round((double)rawTumorAlleleQuality[convert_table[potentialGermlineAllele]]/(double)rawTumorAlleleCount[convert_table[potentialGermlineAllele]]);
					tumorFormat += (to_string(tmp) + ":2");
				}
			}else{
				double homTwotail;
				kt_fisher_exact(coverage, 0, rawNormalAlleleCount[convert_table[potentialGermlineAllele]], rawNormalAlleleCount[convert_table[refLetter]], &homTwotail);
				if(homTwotail>=GERMLINE_VARIANT_FET_CUTOFF){
					std::string tmpS_tmp(1, potentialGermlineAllele);
					altString += ("," + tmpS_tmp);
					
					// normal Format first
					//
					normalFormat = "2/2:" + to_string(rawS[1].length()-rawNormalAlleleCount[convert_table['?']]) + ":" + to_string(rawNormalAlleleCount[convert_table[refLetter]]) + "," + 
											to_string(rawNormalAlleleCount[convert_table[altLetter]]) + "," + to_string(rawNormalAlleleCount[convert_table[potentialGermlineAllele]]) + ":";
					if(rawNormalAlleleCount[convert_table[refLetter]]==0) {
						normalFormat += "0,";
					}
					else {
						int tmp = round((double)rawNormalAlleleQuality[convert_table[refLetter]]/(double)rawNormalAlleleCount[convert_table[refLetter]]);
						normalFormat += (to_string(tmp) + ",");
					}
					if(rawNormalAlleleCount[convert_table[altLetter]]==0) {
						normalFormat += "0,";
					}
					else {
						int tmp = round((double)rawNormalAlleleQuality[convert_table[altLetter]]/(double)rawNormalAlleleCount[convert_table[altLetter]]);
						normalFormat += (to_string(tmp) + ",");
					}
					if(rawNormalAlleleCount[convert_table[potentialGermlineAllele]]==0) {
						normalFormat += "0:.";
					}
					else {
						int tmp = round((double)rawNormalAlleleQuality[convert_table[potentialGermlineAllele]]/(double)rawNormalAlleleCount[convert_table[potentialGermlineAllele]]);
						normalFormat += (to_string(tmp) + ":.");
					}
					// tumor Format
					//
					// always assume the GT is 0/1 unless test result shows otherwise
					//
					// n11       n12         | null (homo-variant)
					// n21       n22         | observation
					//-----------------------+-------------------------------
					// variant   ref         | Total
					//
					double homTwotail;
					int coverage = rawTumorAlleleCount[convert_table[refLetter]] + rawTumorAlleleCount[convert_table[altLetter]];

					kt_fisher_exact(coverage, 0, rawTumorAlleleCount[convert_table[altLetter]], rawTumorAlleleCount[convert_table[refLetter]], &homTwotail);
					
					if(homTwotail>=GERMLINE_VARIANT_FET_CUTOFF) {
						// the test result is not significant, so I cannot reject the hypothesis that this position is homo-variant
						//
						if(rawTumorAlleleCount[convert_table[potentialGermlineAllele]]!=0) {
							tumorFormat = "1/1/2:";
						}
						else {
							tumorFormat = "1/1:";
						}
					}
					else {
						if(rawTumorAlleleCount[convert_table[potentialGermlineAllele]]!=0) {
							tumorFormat = "0/1/2:";
						}
						else {
							tumorFormat = "0/1:";
						}
					}
					tumorFormat += (to_string(rawS[0].length()-rawTumorAlleleCount[convert_table['?']]) + ":" + to_string(rawTumorAlleleCount[convert_table[refLetter]]) + "," + 
									to_string(rawTumorAlleleCount[convert_table[altLetter]]) + "," + to_string(rawTumorAlleleCount[convert_table[potentialGermlineAllele]]) + ":");
					
					if(rawTumorAlleleCount[convert_table[refLetter]]==0) {
						tumorFormat += "0,";
					}
					else {
						int tmp = round((double)rawTumorAlleleQuality[convert_table[refLetter]]/(double)rawTumorAlleleCount[convert_table[refLetter]]);
						tumorFormat += (to_string(tmp) + ",");
					}
					if(rawTumorAlleleCount[convert_table[altLetter]]==0) {
						tumorFormat += "0,";
					}
					else {
						int tmp = round((double)rawTumorAlleleQuality[convert_table[altLetter]]/(double)rawTumorAlleleCount[convert_table[altLetter]]);
						tumorFormat += (to_string(tmp) + ",");
					}
					if(rawTumorAlleleCount[convert_table[potentialGermlineAllele]]==0) {
						tumorFormat += "0:2";
					}
					else {
						int tmp = round((double)rawTumorAlleleQuality[convert_table[potentialGermlineAllele]]/(double)rawTumorAlleleCount[convert_table[potentialGermlineAllele]]);
						tumorFormat += (to_string(tmp) + ":2");
					}
				}else{
					normalFormat = "0/0:" + to_string(rawS[1].length()-rawNormalAlleleCount[convert_table['?']]) + ":" + to_string(rawNormalAlleleCount[convert_table[refLetter]]) + "," + 
											to_string(rawNormalAlleleCount[convert_table[altLetter]]) + ":";
					if(rawNormalAlleleCount[convert_table[refLetter]]==0) {
						normalFormat += "0,";
					}
					else {
						int tmp = round((double)rawNormalAlleleQuality[convert_table[refLetter]]/(double)rawNormalAlleleCount[convert_table[refLetter]]);
						normalFormat += (to_string(tmp) + ",");
					}
					if(rawNormalAlleleCount[convert_table[altLetter]]==0) {
						normalFormat += "0:.";
					}
					else {
						int tmp = round((double)rawNormalAlleleQuality[convert_table[altLetter]]/(double)rawNormalAlleleCount[convert_table[altLetter]]);
						normalFormat += (to_string(tmp) + ":.");
					}
					// tumor Format
					// always assume the GT is 0/1 unless test result shows otherwise
					//
					// n11       n12         | null (homo-variant)
					// n21       n22         | observation
					//-----------------------+-------------------------------
					// variant   ref         | Total
					//
					double homTwotail;
					int coverage = rawTumorAlleleCount[convert_table[refLetter]] + rawTumorAlleleCount[convert_table[altLetter]];

					kt_fisher_exact(coverage, 0, rawTumorAlleleCount[convert_table[altLetter]], rawTumorAlleleCount[convert_table[refLetter]], &homTwotail);
					
					if(homTwotail>=GERMLINE_VARIANT_FET_CUTOFF) {
						// the test result is not significant, so I cannot reject the hypothesis that this position is homo-variant
						//
						tumorFormat = "1/1:";
					}
					else {
						tumorFormat = "0/1:";
					}
					tumorFormat += (to_string(rawS[0].length()-rawTumorAlleleCount[convert_table['?']]) + ":" + to_string(rawTumorAlleleCount[convert_table[refLetter]]) + "," + 
									to_string(rawTumorAlleleCount[convert_table[altLetter]]) + ":");
					
					if(rawTumorAlleleCount[convert_table[refLetter]]==0) {
						tumorFormat += "0,";
					}
					else {
						int tmp = round((double)rawTumorAlleleQuality[convert_table[refLetter]]/(double)rawTumorAlleleCount[convert_table[refLetter]]);
						tumorFormat += (to_string(tmp) + ",");
					}
					if(rawTumorAlleleCount[convert_table[altLetter]]==0) {
						tumorFormat += "0:2";
					}
					else {
						int tmp = round((double)rawTumorAlleleQuality[convert_table[altLetter]]/(double)rawTumorAlleleCount[convert_table[altLetter]]);
						tumorFormat += (to_string(tmp) + ":2");
					}
				}
			}
		}

	}

	output += altString + "\t" + tumorFormat + "\t" + normalFormat + '\n';

}

void WriteHeader(PBLocalFile* outFile, int argc, char * const argv[], bam_hdr_t *tumorH, bam_hdr_t *normalH, const string& refGenome){
	string outf;
	outf += "##MuSE_version=\"" + Version + " Build Date " + buildDate + " Build Time " + buildTime + "\"\n";
	// command line
	//
	outf += "##MuSE_call=\"";
	for(int i = 0; i < argc-1; i++) {
		outf += string(argv[i]) + " ";
	}
	outf += string(argv[argc-1]) + "\"\n";
	// tumor sample
	//
	const char *p = tumorH->text, *tumorSM;
	while((tumorSM = strstr(p, "@RG")) != 0) {
		tumorSM = 0;
		if((tumorSM = strstr(p, "\tSM:")) != 0) tumorSM += 4;
		if(tumorSM) {
			kstring_t kSM;
			memset(&kSM, 0, sizeof(kstring_t));
			char *u;
			int x;
			for(u = (char*)tumorSM; *u && *u != '\t' && *u != '\n'; ++u);
			x = *u; *u = '\0';
			kputs(tumorSM, &kSM);
			if(!kSM.s) {
				std::cerr << "Tumor BAM has no sample information.\n";
				exit(-1);
			}
			outf += string("##TUMOR=\"Sample=") + kSM.s + ",File=" + argv[1] + "\"\n";
			free(kSM.s);
			break;
		} else break;
	}
	// normal sample
	//
	p = normalH->text;
	const char *normalSM;
	while((normalSM = strstr(p, "@RG")) != 0) {
		normalSM = 0;
		if((normalSM = strstr(p, "\tSM:")) != 0) normalSM += 4;
		if(normalSM) {
			kstring_t kSM;
			memset(&kSM, 0, sizeof(kstring_t));
			char *u;
			int x;
			for(u = (char*)normalSM; *u && *u != '\t' && *u != '\n'; ++u);
			x = *u; *u = '\0';
			kputs(normalSM, &kSM);
			if(!kSM.s) {
				std::cerr << "Normal BAM has no sample information.\n";
				exit(-1);
			}
			outf += string("##NORMAL=\"Sample=") + kSM.s + ",File=" + argv[2] + "\"\n";
			free(kSM.s);
			break;
		} else break;
	}
	// contig info
	//
	p = tumorH->text;
	const char *sn, *ln, *as;
	while((sn = strstr(p, "@SQ")) != 0) {
		p = sn + 3;
		sn = ln = as = 0;
		if((sn = strstr(p, "\tSN:")) != 0) sn += 4;
		if((ln = strstr(p, "\tLN:")) != 0) ln += 4;
		if((as = strstr(p, "\tAS:")) != 0) as += 4;
		if(sn && ln) {
			kstring_t kSN, kLN;
			memset(&kSN, 0, sizeof(kstring_t));
			memset(&kLN, 0, sizeof(kstring_t));
			char *u, *v;
			int x, y;
			for(u = (char*)sn; *u && *u != '\t' && *u != '\n'; ++u);
			for(v = (char*)ln; *v && *v != '\t' && *v != '\n'; ++v);
			x = *u; y = *v; *u = *v = '\0';
			kputs(sn, &kSN);
			kputs(ln, &kLN);
			*u = x; *v = y;
			outf += string("##contig=<ID=") + kSN.s + ",length=" + kLN.s;
			free(kSN.s);
			free(kLN.s);
			if(as) {
				kstring_t kAS;
				memset(&kAS, 0, sizeof(kstring_t));
				char *w;
				int z;
				for(w = (char*)as; *w && *w != '\t' && *w != '\n'; ++w);
				z = *w; *w = '\0';
				kputs(as, &kAS);
				*w = z;
				outf += string(",assembly=") + kAS.s + ">\n";
				free(kAS.s);
			}
			else {
				outf += ">\n";
			}
		} else break;
		p = (sn>ln) ? ((sn>as)?sn:as) : ((ln>as)?ln:as);
	}
	// reference file
	//
	outf += "##reference=file://" + refGenome + '\n';
	outFile->pbwrite(outf.c_str(), 1, outf.size());
}

void processPileup(PBReader* reader, mplp_conf_t* conf, std::atomic<bool>& gatherDone, PileUpLockFreeQ& processQ, std::atomic<uint32_t>& processQSize){
	Region* localPileup = nullptr;
	
	const char* ref = nullptr;
	int64_t ref_len = 0;

	WorkMem workMem;
	auto& chrStr = conf->ref->getAllChrStr();

	while(!gatherDone.load() || processQSize.load() > 0){
		if (!processQ.pop(localPileup))
			continue;
		--processQSize;

		auto chrName = reader->mHeader_tumor->target_name[localPileup->tid];
		if (chrStr.find(chrName) != chrStr.end()){
    		ref = conf->ref->getChrStr(chrName).c_str();
			ref_len = conf->ref->getChrLen(chrName);
		}else{
			ref = nullptr;
		}

		if (!ref || localPileup->start >= ref_len){
			localPileup->ifDone.store(true);
			continue;
		}

		workMem.pareparePlp(localPileup->nodes);

		for (int i = 0; i <= min(ref_len, localPileup->end) - localPileup->start; ++i){

			bool ifPush = workMem.getPlp(localPileup->start + i);
			if (ifPush)
				workMem.processPileupFun(reader->mHeader_tumor, conf, localPileup->output, localPileup->tid, i  + localPileup->start, ref, ref_len);
			workMem.clear();

		}

		localPileup->ifDone.store(true);
	}
}