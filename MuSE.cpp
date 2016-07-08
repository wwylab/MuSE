#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <limits.h>
#include <sys/stat.h>
#include <float.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <cctype>
#include <algorithm>
#include <utility>
#include <numeric>
#include <locale>
#include <iterator>

#include "bam.h"
#include "faidx.h"
#include "kstring.h"
#include "sample.h"
#include "bgzf.h"
#include "tabix.h"

#define MPLP_NO_ORPHAN  0x40
#define MPLP_REALN      0x80
#define MPLP_EXT_BAQ    0x800
#define MPLP_ILLUMINA13 0x1000
#define MPLP_IGNORE_RG  0x2000
#define MPLP_PRINT_POS  0x4000
#define MPLP_PRINT_MAPQ 0x8000

#define DEBUG         0
#define PHYML         1
#define POSTERIORMODE 1
#define DREAM         1

#define INPUT_COLUMN_COUNT 30

#define MIN_BRLENS     1.E-8
#define BRLENS_PR_MEAN 0.001

#define MIN_QUALITY_SCORE                5
#define GAP_EVENT_PROXIMITY              5
#define GAP_EVENT_CUTOFF                 3
#define NORMAL_OBSERVATION_COUNT         2
#define NORMAL_OBSERVATION_VAF           0.03
#define NORMAL_OBSERVATION_SUM_QUAL      20
#define GERMLINE_VARIANT_FET_CUTOFF      0.05
#define STRAND_BIAS_FET_CUTOFF           1e-5
#define CLUSTERED_POSITION_MEDIAN_CUTOFF 10
#define CLUSTERED_POSITION_MAD_CUTOFF    3
#define LOW_MAPQ                         10
#define LOW_MEAN_MAPQ                    10
#define HEAVY_SOFTCLIP                   0.3
#define HEAVY_MISMATCH                   100
#define REQUIRED_MAPQ                    30
#define REQUIRED_BASE_QUALITY            25

std::string Version = "v1.0rc";

typedef struct {
	int     flag, min_baseQ, max_depth, depth_cutoff;
	double  minAltFraction, min_output_brlens, normalContRate;
	char    *reg;
	faidx_t *fai;
} mplp_conf_t;

typedef struct {
	bamFile           fp;
	bam_iter_t        iter;
	bam_header_t      *h;
	int               *ref_id;
	int               *ref_len;
	char              **ref_ptr;
	const mplp_conf_t *conf;
} mplp_aux_t;

typedef struct {
	int           n;
	int           *n_plp, *m_plp;
	bam_pileup1_t **plp;
} mplp_pileup_t;

// based on read name, the pair's positions in the pileup
struct pairPosition {
	int first;
	int second;
};

std::vector<int>    weight;
std::vector<double> condLike;
std::vector<double> pMatrix;
std::vector<double> baseFreq;	// if BAYESIAN, prior is Dirichlet(1,1,1,1)
double              kappa;		// if BAYESIAN, prior is Exp(1)
double              brlens;		// if BAYESIAN, prior is Exp(brlenspr_mean)
double              lnlike;
char                refLetter;
char                altLetter;

// debug
unsigned nVisited          = 0;
unsigned nWillBeCalculated = 0;

#if PHYML
std::vector<double> unscaledBaseFreq;
double brlensMin;
double brlensMax;
double min_diff_lk_local;
bool   optimizeBaseFreq;
#endif

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


bool useSecondMutantAllele = false;

// program build information
std::string buildDate = __DATE__;  // e.g. 'Dec 15 2009'
std::string buildTime = __TIME__;  // e.g. '15:25:56'

//================================================================================================= Utilities
std::string IntToString(int x) {
	std::ostringstream o;
    o << x;
    return o.str();
}

int StringToInt(const std::string& s) {
	std::istringstream i(s);
    int x;
    i >> x;
    return x;
}
//================================================================================================= Beta Distribution
#define MY_PI          3.141592653589793238462643383280	/* pi */
#define M_2PI          6.283185307179586476925286766559	/* 2*pi */
#define M_LN_SQRT_2PI  0.918938533204672741780329736406	/* log(sqrt(2*pi)) == log(2*pi)/2 */
#define M_LN_SQRT_PId2 0.225791352644727432363097614947	/* log(sqrt(pi/2)) */
#define M_LOG10_2      0.301029995663981195213738894724	/* log10(2) */
#define MY_LN2         0.693147180559945309417232121458	/* ln(2) */
#define M_1_SQRT_2PI   0.398942280401432677939946059934	/* 1/sqrt(2pi) */
#define M_SQRT_32      5.656854249492380195206754896838	/* sqrt(32) */

#define TRUE       1
#define FALSE      0
#define give_log   log_p
#define ISNAN(x)   isnan(x)
#define FINITE(x)  isfinite(x)
#define ML_NAN     NAN
#define ML_NEGINF -INFINITY
#define ML_POSINF  INFINITY
#define pnorm    pnorm5
#define qnorm    qnorm5
#define dnorm    dnorm4
#define trunc    ftrunc

double r_expm1(double);
double r_log1p(double);

#define repeat                for(;;)
#define RMIN(a,b)             ((a < b)?a:b)
#define RMAX(a,b)             ((a > b)?a:b)
#define R_D__0                (log_p ? ML_NEGINF : 0.)										/* 0 */
#define R_D__1                (log_p ? 0. : 1.)												/* 1 */
#define R_DT_0                (lower_tail ? R_D__0 : R_D__1)									/* 0 */
#define R_DT_1                (lower_tail ? R_D__1 : R_D__0)									/* 1 */
#define R_D_exp(x)            (log_p ? (x) : exp(x))										/* exp(x) */
#define R_DT_qIv(p)           (log_p ? (lower_tail ? exp(p) : - r_expm1(p)) : R_D_Lval(p))	/* #define R_DT_qIv(p)	R_D_Lval(R_D_qIv(p))	*  p  in qF !	*/
#define R_D_val(x)            (log_p ? log(x) : (x))										/*  x  in pF(x,..) */
#define R_D_Lval(p)           (lower_tail ? (p) : (0.5 - (p) + 0.5))						/* Use 0.5 - p + 0.5 to perhaps gain 1 bit of accuracy: p     */
#define R_D_LExp_toms708(x)   (log_p ? R_Log1_Exp_toms708(x) : r_log1p(-x))					/* MM added R_D_LExp, so redefine here in terms of rr_expm1 */
#define R_Log1_Exp_toms708(x) ((x) > -MY_LN2 ? log(-rr_expm1(x)) : r_log1p(-exp(x)))
#define R_D_Lval(p)				(lower_tail ? (p)				: (0.5 - (p) + 0.5)	)			/* Use 0.5 - p + 0.5 to perhaps gain 1 bit of accuracy: p     */
#define R_D_Cval(p)				(lower_tail ? (0.5 - (p) + 0.5) : (p)				)			/* Use 0.5 - p + 0.5 to perhaps gain 1 bit of accuracy: 1 - p */
#define R_DT_CIv(p)           (log_p ? (lower_tail ? -r_expm1(p): exp(p))     : R_D_Cval(p))	/* #define R_DT_CIv(p)	R_D_Cval(R_D_qIv(p))	*  1 - p in qF	*/

#define ME_NONE      0	/*	no error */
#define ME_DOMAIN    1	/*	argument out of domain */
#define ME_RANGE     2	/*	value out of range */
#define ME_NOCONV    4	/*	process did not converge */
#define ME_PRECISION 8	/*	does not have "full" precision */
#define ME_UNDERFLOW 16	/*	and underflow occured (important for IEEE)*/

#define MATHLIB_ERROR(fmt,x)   { printf(fmt,x); exit(1); }
#define MATHLIB_WARNING(fmt,x) printf(fmt,x)
#define ML_ERR_return_NAN      { ML_ERROR(ME_DOMAIN, ""); return ML_NAN; }
#define ML_ERROR(x, s) {										\
if(x > ME_DOMAIN) {												\
	char *msg = "";												\
	switch(x) {													\
		case ME_DOMAIN:											\
			msg = "argument out of domain in '%s'\n";			\
			break;												\
		case ME_RANGE:											\
			msg = "value out of range in '%s'\n";				\
			break;												\
		case ME_NOCONV:											\
			msg = "convergence failed in '%s'\n";				\
			break;												\
		case ME_PRECISION:										\
			msg = "full precision was not achieved in '%s'\n";	\
			break;												\
		case ME_UNDERFLOW:										\
			msg = "underflow occurred in '%s'\n";				\
			break;												\
	}															\
	MATHLIB_WARNING(msg, s);									\
}																\
}

/* Do the boundaries exactly for q*() functions :
 * Often  _LEFT_ = ML_NEGINF , and very often _RIGHT_ = ML_POSINF;
 *
 * R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)  :<==>
 *
 *     R_Q_P01_check(p);
 *     if (p == R_DT_0) return _LEFT_ ;
 *     if (p == R_DT_1) return _RIGHT_;
 *
 * the following implementation should be more efficient (less tests):
 */
#define R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)	\
if (log_p) {									\
	if(p > 0)									\
		return ML_NAN;							\
	if(p == 0)	/* upper bound*/				\
		return lower_tail ? _RIGHT_ : _LEFT_;	\
	if(p == ML_NEGINF)							\
		return lower_tail ? _LEFT_ : _RIGHT_;	\
}												\
else {											\
	if(p < 0 || p > 1)							\
		return ML_NAN;							\
	if(p == 0)									\
		return lower_tail ? _LEFT_ : _RIGHT_;	\
	if(p == 1)									\
		return lower_tail ? _RIGHT_ : _LEFT_;	\
}

typedef enum {
	BUGGY_KINDERMAN_RAMAGE, AHRENS_DIETER, BOX_MULLER, USER_NORM, INVERSION, KINDERMAN_RAMAGE
} N01type;

double gammafn (double);

double fmin2(double x, double y)
{
	if (ISNAN(x)!=0 || ISNAN(y)!=0)
		return x + y;
	return (x < y) ? x : y;
}

double fmax2(double x, double y) {
	if (ISNAN(x)!=0 || ISNAN(y)!=0)
		return x + y;
	return (x < y) ? y : x;
}

double ftrunc(double x) {
	if(x >= 0)
		return floor(x);
	else
		return ceil(x);
}

/*  DESCRIPTION	A version of Marsaglia-MultiCarry */
static unsigned int I1=12345, I2=56789;

void set_seed(unsigned int i1, unsigned int i2)
{
	I1 = i1; I2 = i2;
}
void get_seed(unsigned int *i1, unsigned int *i2)
{
	*i1 = I1; *i2 = I2;
}
double unif_rand(void)
{
	I1= 36969*(I1 & 0177777) + (I1>>16);
	I2= 18000*(I2 & 0177777) + (I2>>16);
	return ((I1 << 16)^(I2 & 0177777)) * 2.328306437080797e-10; /* in [0,1) */
}

/*  DESCRIPTION
 *		Random variates from the standard exponential distribution.
 *  REFERENCE
 *		Ahrens, J.H. and Dieter, U. (1972). Computer methods for sampling from the exponential and normal distributions. Comm. ACM, 15, 873-882.
 */
double exp_rand(void)
{
	/* q[k-1] = sum(log(2)^k / k!)  k=1,..,n, */
	/* The highest n (here 8) is determined by q[n-1] = 1.0 */
	/* within standard precision */
	const static double q[] =
	{
		0.6931471805599453,
		0.9333736875190459,
		0.9888777961838675,
		0.9984959252914960,
		0.9998292811061389,
		0.9999833164100727,
		0.9999985691438767,
		0.9999998906925558,
		0.9999999924734159,
		0.9999999995283275,
		0.9999999999728814,
		0.9999999999985598,
		0.9999999999999289,
		0.9999999999999968,
		0.9999999999999999,
		1.0000000000000000
	};
	double a, u, ustar, umin;
	int i;
	
	a = 0.;
	/* precaution if u = 0 is ever returned */
	u = unif_rand();
	while(u <= 0.0 || u >= 1.0)
		u = unif_rand();
	for (;;) {
		u += u;
		if (u > 1.0)
			break;
		a += q[0];
	}
	u -= 1.;
	
	if (u <= q[0])
		return a + u;
	
	i = 0;
	ustar = unif_rand();
	umin = ustar;
	do {
		ustar = unif_rand();
		if (ustar < umin)
			umin = ustar;
		i++;
	} while (u > q[i]);
	return a + umin * q[0];
}

/*  DESCRIPTION
 *		Compute the density of the normal distribution.
 */
double dnorm4(double x, double mu, double sigma, int give_log)
{
	if (ISNAN(x)!=0 || ISNAN(mu)!=0 || ISNAN(sigma)!=0)
		return x + mu + sigma;
	if(FINITE(sigma)==0)
		return R_D__0;
	if(FINITE(x)==0 && mu == x)
		return ML_NAN;/* x-mu is NaN */
	if (sigma <= 0) {
		if (sigma < 0)
			ML_ERR_return_NAN
		/* sigma == 0 */
			return (x == mu) ? ML_POSINF : R_D__0;
	}
	x = (x - mu) / sigma;
	
	if(FINITE(x)==0)
		return R_D__0;
	return (give_log ? -(M_LN_SQRT_2PI + 0.5*x*x + log(sigma)) : M_1_SQRT_2PI * exp(-0.5*x*x) / sigma);
}

/*  DESCRIPTION
 *		The main computation evaluates near-minimax approximations derived from those in "Rational Chebyshev approximations for the error
 *		function" by W. J. Cody, Math. Comp., 1969, 631-637.  This transportable program uses rational functions that theoretically
 *		approximate the normal distribution function to at least 18 significant decimal digits.  The accuracy achieved depends on the
 *		arithmetic system, the compiler, the intrinsic functions, and proper selection of the machine-dependent constants.
 *  REFERENCE
 *		Cody, W. D. (1993). ALGORITHM 715: SPECFUN - A Portable FORTRAN Package of Special Function Routines and Test Drivers". ACM Transactions on Mathematical Software. 19, 22-32.
 *  EXTENSIONS
 *		The "_both" , lower, upper, and log_p  variants were added by Martin Maechler, Jan.2000; as well as r_log1p() and similar improvements later on.
 *
 *		James M. Rath contributed bug report PR#699 and patches correcting SIXTEN and if() clauses {with a bug: "|| instead of &&" -> PR #2883) more in line with the original Cody code.
 */
#define SIXTEN	16 /* Cutoff allowing exact "*" and "/" */

void pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p)
{
	/* i_tail in {0,1,2} means: "lower", "upper", or "both" :
	 if(lower) return  *cum := P[X <= x]
	 if(upper) return *ccum := P[X >  x] = 1 - P[X <= x]
	 */
	const static double a[5] = {
		2.2352520354606839287,
		161.02823106855587881,
		1067.6894854603709582,
		18154.981253343561249,
		0.065682337918207449113
	};
	const static double b[4] = {
		47.20258190468824187,
		976.09855173777669322,
		10260.932208618978205,
		45507.789335026729956
	};
	const static double c[9] = {
		0.39894151208813466764,
		8.8831497943883759412,
		93.506656132177855979,
		597.27027639480026226,
		2494.5375852903726711,
		6848.1904505362823326,
		11602.651437647350124,
		9842.7148383839780218,
		1.0765576773720192317e-8
	};
	const static double d[8] = {
		22.266688044328115691,
		235.38790178262499861,
		1519.377599407554805,
		6485.558298266760755,
		18615.571640885098091,
		34900.952721145977266,
		38912.003286093271411,
		19685.429676859990727
	};
	const static double p[6] = {
		0.21589853405795699,
		0.1274011611602473639,
		0.022235277870649807,
		0.001421619193227893466,
		2.9112874951168792e-5,
		0.02307344176494017303
	};
	const static double q[5] = {
		1.28426009614491121,
		0.468238212480865118,
		0.0659881378689285515,
		0.00378239633202758244,
		7.29751555083966205e-5
	};
	
	double xden, xnum, temp, del, eps, xsq, y;
#ifdef NO_DENORMS
	double min = DBL_MIN;
#endif
	int i, lower, upper;
	
	if(ISNAN(x)!=0) { *cum = *ccum = x; return; }
	
	/* Consider changing these : */
	eps = DBL_EPSILON * 0.5;
	
	/* i_tail in {0,1,2} =^= {lower, upper, both} */
	lower = i_tail != 1;
	upper = i_tail != 0;
	
	y = fabs(x);
	if (y <= 0.67448975) { /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
		if (y > eps) {
			xsq  = x * x;
			xnum = a[4] * xsq;
			xden = xsq;
			for (i = 0; i < 3; ++i) {
				xnum = (xnum + a[i]) * xsq;
				xden = (xden + b[i]) * xsq;
			}
		} else xnum = xden = 0.0;
		
		temp = x * (xnum + a[3]) / (xden + b[3]);
		if(lower)  *cum  = 0.5 + temp;
		if(upper)  *ccum = 0.5 - temp;
		if(log_p) {
			if(lower)  *cum  = log(*cum);
			if(upper)  *ccum = log(*ccum);
		}
	}
	else if (y <= M_SQRT_32) {
		
		/* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */
		
		xnum = c[8] * y;
		xden = y;
		for (i = 0; i < 7; ++i) {
			xnum = (xnum + c[i]) * y;
			xden = (xden + d[i]) * y;
		}
		temp = (xnum + c[7]) / (xden + d[7]);
		
		xsq = trunc(y * SIXTEN) / SIXTEN;
		del = (y - xsq) * (y + xsq);
		if(log_p) {
			*cum = (-xsq * xsq * 0.5) + (-del * 0.5) + log(temp);
			if((lower && x > 0.) || (upper && x <= 0.))
				*ccum = r_log1p(-exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * temp);
		}
		else {
			*cum  = exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * temp;
			*ccum = 1.0 - *cum;
		}
		
		if (x > 0.) {       // swap  ccum <--> cum
			temp = *cum;
			if(lower)
				*cum = *ccum;
			*ccum = temp;
		}
	}
	
	/* else	  |x| > sqrt(32) = 5.657 :
	 * the next two case differentiations were really for lower=T, log=F
	 * Particularly	 *not*	for  log_p !
	 
	 * Cody had (-37.5193 < x  &&  x < 8.2924) ; R originally had y < 50
	 *
	 * Note that we do want symmetry(0), lower/upper -> hence use y
	 */
	else if(log_p
			/*  ^^^^^ MM FIXME: can speedup for log_p and much larger |x| !
			 * Then, make use of  Abramowitz & Stegun, 26.2.13, something like
			 
			 xsq = x*x;
			 
			 if(xsq * DBL_EPSILON < 1.)
			 del = (1. - (1. - 5./(xsq+6.)) / (xsq+4.)) / (xsq+2.);
			 else
			 del = 0.;
			 *cum = -.5*xsq - M_LN_SQRT_2PI - log(x) + r_log1p(-del);
			 *ccum = r_log1p(-exp(*cum)); /.* ~ log(1) = 0 *./
			 
			 swap_tail;
			 
			 */
			|| (lower && -37.5193 < x  &&  x < 8.2924)
			|| (upper && -8.2924  < x  &&  x < 37.5193)
			) {
		
		/* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
		xsq = 1.0 / (x * x);
		xnum = p[5] * xsq;
		xden = xsq;
		for (i = 0; i < 4; ++i) {
			xnum = (xnum + p[i]) * xsq;
			xden = (xden + q[i]) * xsq;
		}
		temp = xsq * (xnum + p[4]) / (xden + q[4]);
		temp = (M_1_SQRT_2PI - temp) / y;
		
		xsq = trunc(x * SIXTEN) / SIXTEN;
		del = (x - xsq) * (x + xsq);
		if(log_p) {
			*cum = (-xsq * xsq * 0.5) + (-del * 0.5) + log(temp);
			if((lower && x > 0.) || (upper && x <= 0.))
				*ccum = r_log1p(-exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * temp);
		}
		else {
			*cum  = exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * temp;
			*ccum = 1.0 - *cum;
		}
		if (x > 0.) {      // swap  ccum <--> cum
			temp = *cum;
			if(lower)
				*cum = *ccum;
			*ccum = temp;
		}
	}
	else { /* no log_p , large x such that probs are 0 or 1 */
		if(x > 0) {
			*cum = 1.;
			*ccum = 0.;
		}
		else {
			*cum = 0.;
			*ccum = 1.;
		}
	}
	
#ifdef NO_DENORMS
	/* do not return "denormalized" -- we do in R */
	if(log_p) {
		if(*cum > -min)	 *cum = -0.;
		if(*ccum > -min)*ccum = -0.;
	}
	else {
		if(*cum < min)	 *cum = 0.;
		if(*ccum < min)	*ccum = 0.;
	}
#endif
	
	return;
}

double pnorm5(double x, double mu, double sigma, int lower_tail, int log_p)
{
	double p, cp;
	
	/* Note: The structure of these checks has been carefully thought through.
	 * For example, if x == mu and sigma == 0, we get the correct answer 1.
	 */
	if(ISNAN(x)!=0 || ISNAN(mu)!=0 || ISNAN(sigma)!=0)
		return x + mu + sigma;
	if(FINITE(x)==0 && mu == x)
		return ML_NAN;/* x-mu is NaN */
	if (sigma <= 0) {
		if(sigma < 0)
			ML_ERR_return_NAN
		/* sigma = 0 : */
			return (x < mu) ? R_DT_0 : R_DT_1;
	}
	p = (x - mu) / sigma;
	if(FINITE(p)==0)
		return (x < mu) ? R_DT_0 : R_DT_1;
	x = p;
	
	pnorm_both(x, &p, &cp, (lower_tail ? 0 : 1), log_p);
	
	return(lower_tail ? p : cp);
}

/*  DESCRIPTION
 *		Compute the quantile function for the normal distribution.
 *		For small to moderate probabilities, algorithm referenced below is used to obtain an initial approximation which is polished with a final Newton step.
 *		For very large arguments, an algorithm of Wichura is used.
 *  REFERENCE
 *		Beasley, J. D. and S. G. Springer (1977). Algorithm AS 111: The percentage points of the normal distribution. Applied Statistics, 26, 118-121.
 *		Wichura, M.J. (1988). Algorithm AS 241: The Percentage Points of the Normal Distribution. Applied Statistics, 37, 477-484.
 */
double qnorm5(double p, double mu, double sigma, int lower_tail, int log_p)
{
	double p_, q, r, val;
	
	if (ISNAN(p)!=0 || ISNAN(mu)!=0 || ISNAN(sigma)!=0)
		return p + mu + sigma;
	
	R_Q_P01_boundaries(p, ML_NEGINF, ML_POSINF);
	
	if(sigma  < 0)	ML_ERR_return_NAN;
	if(sigma == 0)	return mu;
	
	p_ = R_DT_qIv(p);/* real lower_tail prob. p */
	q  = p_ - 0.5;
	
	/*-- use AS 241 --- */
	/* double ppnd16_(double *p, long *ifault)*/
	/*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
	 *      Produces the normal deviate Z corresponding to a given lower tail area of P; Z is accurate to about 1 part in 10**16.
	 *      (original fortran code used PARAMETER(..) for the coefficients and provided hash codes for checking them...)
	 */
	if (fabs(q) <= .425) {/* 0.075 <= p <= 0.925 */
		r = .180625 - q * q;
		val =
		q * (((((((r * 2509.0809287301226727 +
				   33430.575583588128105) * r + 67265.770927008700853) * r +
				 45921.953931549871457) * r + 13731.693765509461125) * r +
			   1971.5909503065514427) * r + 133.14166789178437745) * r +
			 3.387132872796366608)
		/ (((((((r * 5226.495278852854561 +
				 28729.085735721942674) * r + 39307.89580009271061) * r +
			   21213.794301586595867) * r + 5394.1960214247511077) * r +
			 687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
	}
	else { /* closer than 0.075 from {0,1} boundary */
		
		/* r = min(p, 1-p) < 0.075 */
		if (q > 0)
			r = R_DT_CIv(p);/* 1-p */
		else
			r = p_;/* = R_DT_Iv(p) ^=  p */
		
		r = sqrt(- ((log_p && ((lower_tail && q <= 0) || (!lower_tail && q > 0))) ? p : /* else */ log(r))); /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */
		
		if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
			r += -1.6;
			val = (((((((r * 7.7454501427834140764e-4 +
						 .0227238449892691845833) * r + .24178072517745061177) *
					   r + 1.27045825245236838258) * r +
					  3.64784832476320460504) * r + 5.7694972214606914055) *
					r + 4.6303378461565452959) * r +
				   1.42343711074968357734)
			/ (((((((r *
					 1.05075007164441684324e-9 + 5.475938084995344946e-4) *
					r + .0151986665636164571966) * r +
				   .14810397642748007459) * r + .68976733498510000455) *
				 r + 1.6763848301838038494) * r +
				2.05319162663775882187) * r + 1.);
		}
		else { /* very close to  0 or 1 */
			r += -5.;
			val = (((((((r * 2.01033439929228813265e-7 +
						 2.71155556874348757815e-5) * r +
						.0012426609473880784386) * r + .026532189526576123093) *
					  r + .29656057182850489123) * r +
					 1.7848265399172913358) * r + 5.4637849111641143699) *
				   r + 6.6579046435011037772)
			/ (((((((r *
					 2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
					r + 1.8463183175100546818e-5) * r +
				   7.868691311456132591e-4) * r + .0148753612908506148525)
				 * r + .13692988092273580531) * r +
				.59983220655588793769) * r + 1.);
		}
		
		if(q < 0.0)
			val = -val;
		/* return (q >= 0.)? r : -r ;*/
	}
	return mu + sigma * val;
}

/*  DESCRIPTION
 *		Random variates from the STANDARD normal distribution  N(0,1). Is called from  rnorm(..), but also rt(), rf(), rgamma(), ...
 *  REFERENCE
 *    Ahrens, J.H. and Dieter, U. Extensions of Forsythe's method for random sampling from the normal distribution. Math. Comput. 27, 927-937.
 *	NOTES
 *    The definitions of the constants a[k], d[k], t[k] and h[k] are according to the abovementioned article
 */
static double BM_norm_keep = 0.0;
N01type N01_kind = INVERSION; /* Different kinds of "N(0,1)" generators :*/

double norm_rand(void)
{
	const static double a[32] =
	{
		0.0000000, 0.03917609, 0.07841241, 0.1177699,
		0.1573107, 0.19709910, 0.23720210, 0.2776904,
		0.3186394, 0.36012990, 0.40225010, 0.4450965,
		0.4887764, 0.53340970, 0.57913220, 0.6260990,
		0.6744898, 0.72451440, 0.77642180, 0.8305109,
		0.8871466, 0.94678180, 1.00999000, 1.0775160,
		1.1503490, 1.22985900, 1.31801100, 1.4177970,
		1.5341210, 1.67594000, 1.86273200, 2.1538750
	};
	
	const static double d[31] =
	{
		0.0000000, 0.0000000, 0.0000000, 0.0000000,
		0.0000000, 0.2636843, 0.2425085, 0.2255674,
		0.2116342, 0.1999243, 0.1899108, 0.1812252,
		0.1736014, 0.1668419, 0.1607967, 0.1553497,
		0.1504094, 0.1459026, 0.1417700, 0.1379632,
		0.1344418, 0.1311722, 0.1281260, 0.1252791,
		0.1226109, 0.1201036, 0.1177417, 0.1155119,
		0.1134023, 0.1114027, 0.1095039
	};
	
	const static double t[31] =
	{
		7.673828e-4, 0.002306870, 0.003860618, 0.005438454,
		0.007050699, 0.008708396, 0.010423570, 0.012209530,
		0.014081250, 0.016055790, 0.018152900, 0.020395730,
		0.022811770, 0.025434070, 0.028302960, 0.031468220,
		0.034992330, 0.038954830, 0.043458780, 0.048640350,
		0.054683340, 0.061842220, 0.070479830, 0.081131950,
		0.094624440, 0.112300100, 0.136498000, 0.171688600,
		0.227624100, 0.330498000, 0.584703100
	};
	
	const static double h[31] =
	{
		0.03920617, 0.03932705, 0.03950999, 0.03975703,
		0.04007093, 0.04045533, 0.04091481, 0.04145507,
		0.04208311, 0.04280748, 0.04363863, 0.04458932,
		0.04567523, 0.04691571, 0.04833487, 0.04996298,
		0.05183859, 0.05401138, 0.05654656, 0.05953130,
		0.06308489, 0.06737503, 0.07264544, 0.07926471,
		0.08781922, 0.09930398, 0.11555990, 0.14043440,
		0.18361420, 0.27900160, 0.70104740
	};
	
	/*----------- Constants and definitions for  Kinderman - Ramage --- */
	/*
	 *  REFERENCE
	 *		Kinderman A. J. and Ramage J. G. (1976). Computer generation of normal random variables. JASA 71, 893-896.
	 */
	
#define C1		0.398942280401433
#define C2		0.180025191068563
#define g(x)	(C1*exp(-x*x/2.0)-C2*(A-x))
	
	const static double A = 2.216035867166471;
	
	double s, u1, w, y, u2, u3, aa, tt, theta, R;
	int i;
	
	switch(N01_kind) {
			
		case  AHRENS_DIETER: /* see Reference above */
			u1 = unif_rand();
			s = 0.0;
			if (u1 > 0.5)
				s = 1.0;
			u1 = u1 + u1 - s;
			u1 *= 32.0;
			i = (int) u1;
			if (i == 32)
				i = 31;
			if (i != 0) {
				u2 = u1 - i;
				aa = a[i - 1];
				while (u2 <= t[i - 1]) {
					u1 = unif_rand();
					w  = u1 * (a[i] - aa);
					tt = (w * 0.5 + aa) * w;
					repeat {
						if (u2 > tt)
							goto deliver;
						u1 = unif_rand();
						if (u2 < u1)
							break;
						tt = u1;
						u2 = unif_rand();
					}
					u2 = unif_rand();
				}
				w = (u2 - t[i - 1]) * h[i - 1];
			}
			else {
				i = 6;
				aa = a[31];
				repeat {
					u1 = u1 + u1;
					if (u1 >= 1.0)
						break;
					aa = aa + d[i - 1];
					i = i + 1;
				}
				u1 = u1 - 1.0;
				repeat {
					w = u1 * d[i - 1];
					tt = (w * 0.5 + aa) * w;
					repeat {
						u2 = unif_rand();
						if (u2 > tt)
							goto jump;
						u1 = unif_rand();
						if (u2 < u1)
							break;
						tt = u1;
					}
					u1 = unif_rand();
				}
			jump:;
			}
			
		deliver:
			y = aa + w;
			return (s == 1.0) ? -y : y;
			
		case BUGGY_KINDERMAN_RAMAGE: /* see Reference above */
			/* NOTE: this has problems, but is retained for reproducibility of older codes, with the same numeric code */
			u1 = unif_rand();
			if(u1 < 0.884070402298758) {
				u2 = unif_rand();
				return A*(1.13113163544180*u1+u2-1);
			}
			
			if(u1 >= 0.973310954173898) { /* tail: */
				repeat {
					u2 = unif_rand();
					u3 = unif_rand();
					tt = (A*A-2*log(u3));
					if( u2*u2<(A*A)/tt )
						return (u1 < 0.986655477086949) ? sqrt(tt) : -sqrt(tt);
				}
			}
			
			if(u1 >= 0.958720824790463) { /* region3: */
				repeat {
					u2 = unif_rand();
					u3 = unif_rand();
					tt = A - 0.630834801921960* fmin2(u2,u3);
					if(fmax2(u2,u3) <= 0.755591531667601)
						return (u2<u3) ? tt : -tt;
					if(0.034240503750111*fabs(u2-u3) <= g(tt))
						return (u2<u3) ? tt : -tt;
				}
			}
			
			if(u1 >= 0.911312780288703) { /* region2: */
				repeat {
					u2 = unif_rand();
					u3 = unif_rand();
					tt = 0.479727404222441+1.105473661022070*fmin2(u2,u3);
					if( fmax2(u2,u3)<=0.872834976671790 )
						return (u2<u3) ? tt : -tt;
					if( 0.049264496373128*fabs(u2-u3)<=g(tt) )
						return (u2<u3) ? tt : -tt;
				}
			}
			
			/* ELSE	 region1: */
			repeat {
				u2 = unif_rand();
				u3 = unif_rand();
				tt = 0.479727404222441-0.595507138015940*fmin2(u2,u3);
				if(fmax2(u2,u3) <= 0.805577924423817)
					return (u2<u3) ? tt : -tt;
			}
			
		case BOX_MULLER:
			if(BM_norm_keep != 0.0) { /* An exact test is intentional */
				s = BM_norm_keep;
				BM_norm_keep = 0.0;
				return s;
			} else {
				theta = 2 * MY_PI * unif_rand();
				R = sqrt(-2 * log(unif_rand())) + 10*DBL_MIN; /* ensure non-zero */
				BM_norm_keep = R * sin(theta);
				return R * cos(theta);
			}
			//#ifndef MATHLIB_STANDALONE
			//    case USER_NORM:
			//	return *((double *) User_norm_fun());
			//#endif
		case INVERSION:
#define BIG 134217728 /* 2^27 */
			/* unif_rand() alone is not of high enough precision */
			u1 = unif_rand();
			u1 = (int)(BIG*u1) + unif_rand();
			return qnorm5(u1/BIG, 0.0, 1.0, 1, 0);
			
		case KINDERMAN_RAMAGE: /* see Reference above */
			/* corrected version from Josef Leydold */
			u1 = unif_rand();
			if(u1 < 0.884070402298758) {
				u2 = unif_rand();
				return A*(1.131131635444180*u1+u2-1);
			}
			
			if(u1 >= 0.973310954173898) { /* tail: */
				repeat {
					u2 = unif_rand();
					u3 = unif_rand();
					tt = (A*A-2*log(u3));
					if( u2*u2<(A*A)/tt )
						return (u1 < 0.986655477086949) ? sqrt(tt) : -sqrt(tt);
				}
			}
			
			if(u1 >= 0.958720824790463) { /* region3: */
				repeat {
					u2 = unif_rand();
					u3 = unif_rand();
					tt = A - 0.630834801921960* fmin2(u2,u3);
					if(fmax2(u2,u3) <= 0.755591531667601)
						return (u2<u3) ? tt : -tt;
					if(0.034240503750111*fabs(u2-u3) <= g(tt))
						return (u2<u3) ? tt : -tt;
				}
			}
			
			if(u1 >= 0.911312780288703) { /* region2: */
				repeat {
					u2 = unif_rand();
					u3 = unif_rand();
					tt = 0.479727404222441+1.105473661022070*fmin2(u2,u3);
					if( fmax2(u2,u3)<=0.872834976671790 )
						return (u2<u3) ? tt : -tt;
					if( 0.049264496373128*fabs(u2-u3)<=g(tt) )
						return (u2<u3) ? tt : -tt;
				}
			}
			/* ELSE	 region1: */
			repeat {
				u2 = unif_rand();
				u3 = unif_rand();
				tt = 0.479727404222441-0.595507138015940*fmin2(u2,u3);
				if (tt < 0.)
					continue;
				if(fmax2(u2,u3) <= 0.805577924423817)
					return (u2<u3) ? tt : -tt;
				if(0.053377549506886*fabs(u2-u3) <= g(tt))
					return (u2<u3) ? tt : -tt;
			}
			
		default:
			MATHLIB_ERROR("norm_rand(): invalid N01_kind: %d\n", N01_kind)
			return 0.0;/*- -Wall */
	}/*switch*/
}

/*  DESCRIPTION
 *		Evaluates the "deviance part"
 *			bd0(x,M) :=  M * D0(x/M) = M*[ x/M * log(x/M) + 1 - (x/M) ] =  x * log(x/M) + M - x
 *		where M = E[X] = n*p (or = lambda), for	  x, M > 0
 *		in a manner that should be stable (with small relative error) for all x and M=np. In particular for x/np close to 1, direct
 *		evaluation fails, and evaluation is based on the Taylor series of log((1+v)/(1-v)) with v = (x-np)/(x+np).
 */
double bd0(double x, double np) {
	double ej, s, s1, v;
	int j;
	
	if(FINITE(x)==0 || FINITE(np)==0 || np == 0.0)
		ML_ERR_return_NAN
		
		if (fabs(x-np) < 0.1*(x+np)) {
			v = (x-np)/(x+np);
			s = (x-np)*v;/* s using v -- change by MM */
			ej = 2*x*v;
			v = v*v;
			for (j=1; ; j++) { /* Taylor series */
				ej *= v;
				s1 = s+ej/((j<<1)+1);
				if (s1==s) /* last term was effectively 0 */
					return(s1);
				s = s1;
			}
		}
	/* else:  | x - np |  is not too small */
	return(x*log(x/np)+np-x);
}

/*  DESCRIPTION
 *		"chebyshev_init" determines the number of terms for the double precision orthogonal series "dos" needed to insure the error is no larger than "eta".
 *		Ordinarily eta will be chosen to be one-tenth machine precision.
 *		"chebyshev_eval" evaluates the n-term Chebyshev series "a" at "x".
 *	NOTES
 *		These routines are translations into C of Fortran routines by W. Fullerton of Los Alamos Scientific Laboratory.
 *		Based on the Fortran routine dcsevl by W. Fullerton. Adapted from R. Broucke, Algorithm 446, CACM., 16, 254 (1973).
 */
double chebyshev_eval(double x, const double *a, const int n) {
	double b0, b1, b2, twox;
	int i;
	
	if (n < 1 || n > 1000)
		ML_ERR_return_NAN
		
		if (x < -1.1 || x > 1.1)
			ML_ERR_return_NAN
			
			twox = x * 2;
	b2 = b1 = 0;
	b0 = 0;
	for (i = 1; i <= n; i++) {
		b2 = b1;
		b1 = b0;
		b0 = twox * b1 - b2 + a[n - i];
	}
	return (b0 - b2) * 0.5;
}

/*  DESCRIPTION
 *		Compute the relative error logarithm.
 *			log(1 + x)
 *	NOTES
 *		This code is a translation of the Fortran subroutine `dlnrel' written by W. Fullerton of Los Alamos Scientific Laboratory.
 */
double r_log1p(double x) {
	/* series for r_log1p on the interval -.375 to .375
	 *				     with weighted error   6.35e-32
	 *				      log weighted error  31.20
	 *			    significant figures required  30.93
	 *				 decimal places required  32.01
	 */
	const static double alnrcs[43] = {
		+.10378693562743769800686267719098e+1,
		-.13364301504908918098766041553133e+0,
		+.19408249135520563357926199374750e-1,
		-.30107551127535777690376537776592e-2,
		+.48694614797154850090456366509137e-3,
		-.81054881893175356066809943008622e-4,
		+.13778847799559524782938251496059e-4,
		-.23802210894358970251369992914935e-5,
		+.41640416213865183476391859901989e-6,
		-.73595828378075994984266837031998e-7,
		+.13117611876241674949152294345011e-7,
		-.23546709317742425136696092330175e-8,
		+.42522773276034997775638052962567e-9,
		-.77190894134840796826108107493300e-10,
		+.14075746481359069909215356472191e-10,
		-.25769072058024680627537078627584e-11,
		+.47342406666294421849154395005938e-12,
		-.87249012674742641745301263292675e-13,
		+.16124614902740551465739833119115e-13,
		-.29875652015665773006710792416815e-14,
		+.55480701209082887983041321697279e-15,
		-.10324619158271569595141333961932e-15,
		+.19250239203049851177878503244868e-16,
		-.35955073465265150011189707844266e-17,
		+.67264542537876857892194574226773e-18,
		-.12602624168735219252082425637546e-18,
		+.23644884408606210044916158955519e-19,
		-.44419377050807936898878389179733e-20,
		+.83546594464034259016241293994666e-21,
		-.15731559416479562574899253521066e-21,
		+.29653128740247422686154369706666e-22,
		-.55949583481815947292156013226666e-23,
		+.10566354268835681048187284138666e-23,
		-.19972483680670204548314999466666e-24,
		+.37782977818839361421049855999999e-25,
		-.71531586889081740345038165333333e-26,
		+.13552488463674213646502024533333e-26,
		-.25694673048487567430079829333333e-27,
		+.48747756066216949076459519999999e-28,
		-.92542112530849715321132373333333e-29,
		+.17578597841760239233269760000000e-29,
		-.33410026677731010351377066666666e-30,
		+.63533936180236187354180266666666e-31,
	};
	
#ifdef NOMORE_FOR_THREADS
	static int nlnrel = 0;
	static double xmin = 0.0;
	
	if (xmin == 0.0)
		xmin = -1 + sqrt(DBL_EPSILON);/*was sqrt(d1mach(4)); */
	if (nlnrel == 0) /* initialize chebychev coefficients */
		nlnrel = chebyshev_init(alnrcs, 43, DBL_EPSILON/20);/*was .1*d1mach(3)*/
#else
# define nlnrel 22
	const static double xmin = -0.999999985;
	/* 22: for IEEE double precision where DBL_EPSILON =  2.22044604925031e-16 */
#endif
	
	if (x == 0.) return 0.;/* speed */
	if (x == -1) return(ML_NEGINF);
	if (x  < -1) ML_ERR_return_NAN
		
		if (fabs(x) <= .375) {
			/* Improve on speed (only); again give result accurate to IEEE double precision: */
			if(fabs(x) < .5 * DBL_EPSILON)
				return x;
			
			if( (0 < x && x < 1e-8) || (-1e-9 < x && x < 0))
				return x * (1 - .5 * x);
			/* else */
			return x * (1 - x * chebyshev_eval(x / .375, alnrcs, nlnrel));
		}
	/* else */
	if (x < xmin) {
		/* answer less than half precision because x too near -1 */
		ML_ERROR(ME_PRECISION, "r_log1p");
	}
	return log(1 + x);
}

static double rlog1(double x) {
	/* -----------------------------------------------------------------------
	 *             Evaluation of the function  x - ln(1 + x)
	 * ----------------------------------------------------------------------- */
	
	static double a  = .0566749439387324;
	static double b  = .0456512608815524;
	static double p0 = .333333333333333;
	static double p1 = -.224696413112536;
	static double p2 = .00620886815375787;
	static double q1 = -1.27408923933623;
	static double q2 = .354508718369557;
	
	double h, r, t, w, w1;
	
	if (x < -0.39 || x > 0.57) { /* direct evaluation */
		w = x + 0.5 + 0.5;
		return x - log(w);
	}
	/* else */
	if (x < -0.18) { /* L10: */
		h = x + .3;
		h /= .7;
		w1 = a - h * .3;
	}
	else if (x > 0.18) { /* L20: */
		h = x * .75 - .25;
		w1 = b + h / 3.0;
	}
	else { /*		Argument Reduction */
		h = x;
		w1 = 0.0;
	}
	
	/* L30:              	Series Expansion */
	
	r = h / (h + 2.0);
	t = r * r;
	w = ((p2 * t + p1) * t + p0) / ((q2 * t + q1) * t + 1.0);
	return t * 2.0 * (1.0 / (1.0 - r) - r * w) + w1;
}

/*  DESCRIPTION
 *		Compute the Exponential minus 1
 *			exp(x) - 1
 *		accurately also when x is close to zero, i.e. |x| << 1
 *  NOTES
 *	As r_log1p(), this is a standard function in some C libraries, particularly GNU and BSD (but is neither ISO/ANSI C nor POSIX).
 *	We supply a substitute for the case when there is no system one.
 */
double r_expm1(double x) {
	double y, a = fabs(x);
	if (a < DBL_EPSILON)	return x;
	if (a > 0.697)			return exp(x) - 1;  /* negligible cancellation */
	if (a > 1e-8)
		y = exp(x) - 1;
	else /* Taylor expansion, more accurate in this range */
		y = (x / 2 + 1) * x;
	
	/* Newton step for solving   log(1 + y) = x   for y : */
	/* WARNING: does not work for y ~ -1: bug in 1.5.0 */
	y -= (1 + y) * (r_log1p (y) - x);
	return y;
}

double rr_expm1(double x) {
	/* ----------------------------------------------------------------------- */
	/*            EVALUATION OF THE FUNCTION EXP(X) - 1 */
	/* ----------------------------------------------------------------------- */
	
	static double p1 = 9.14041914819518e-10;
	static double p2 = .0238082361044469;
	static double q1 = -.499999999085958;
	static double q2 = .107141568980644;
	static double q3 = -.0119041179760821;
	static double q4 = 5.95130811860248e-4;
	
	if (fabs(x) <= 0.15) {
		return x * (((p2 * x + p1) * x + 1.0) / ((((q4 * x + q3) * x + q2) * x + q1) * x + 1.0));
	}
	else { /* |x| > 0.15 : */
		double w = exp(x);
		if (x > 0.0)
			return w * (0.5 - 1.0 / w + 0.5);
		else
			return w - 0.5 - 0.5;
	}
}

double Rf_d1mach(int i) {
	switch(i) {
		case 1:
			return DBL_MIN;
		case 2:
			return DBL_MAX;
		case 3: /* = FLT_RADIX  ^ - DBL_MANT_DIG for IEEE:  = 2^-53 = 1.110223e-16 = .5*DBL_EPSILON */
			return 0.5*DBL_EPSILON;
		case 4: /* = FLT_RADIX  ^ (1- DBL_MANT_DIG) = for IEEE:  = 2^-52 = DBL_EPSILON */
			return DBL_EPSILON;
		case 5:
			return M_LOG10_2;
		default:
			return 0.0;
	}
}

int Rf_i1mach(int i) {
	switch(i) {
		case  1: return 5;
		case  2: return 6;
		case  3: return 0;
		case  4: return 0;
		case  5: return CHAR_BIT * sizeof(int);
		case  6: return sizeof(int)/sizeof(char);
		case  7: return 2;
		case  8: return CHAR_BIT * sizeof(int) - 1;
		case  9: return INT_MAX;
		case 10: return FLT_RADIX;
		case 11: return FLT_MANT_DIG;
		case 12: return FLT_MIN_EXP;
		case 13: return FLT_MAX_EXP;
		case 14: return DBL_MANT_DIG;
		case 15: return DBL_MIN_EXP;
		case 16: return DBL_MAX_EXP;
		default: return 0;
	}
}

static double exparg(int l) {
	/* -------------------------------------------------------------------- */
	/*     IF L = 0 THEN  EXPARG(L) = THE LARGEST POSITIVE W FOR WHICH */
	/*     EXP(W) CAN BE COMPUTED. */
	
	/*     IF L IS NONZERO THEN  EXPARG(L) = THE LARGEST NEGATIVE W FOR */
	/*     WHICH THE COMPUTED VALUE OF EXP(W) IS NONZERO. */
	
	/*     NOTE... ONLY AN APPROXIMATE VALUE FOR EXPARG(L) IS NEEDED. */
	/* -------------------------------------------------------------------- */
	
	static double const lnb = .69314718055995;
	int m;
	
	if (l == 0) {
		m = Rf_i1mach(16);
		return m * lnb * .99999;
	}
	m = Rf_i1mach(15) - 1;
	return m * lnb * .99999;
}

static double erfc1(int ind, double x) {
	/* ----------------------------------------------------------------------- */
	/*         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION */
	
	/*          ERFC1(IND,X) = ERFC(X)            IF IND = 0 */
	/*          ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE */
	/* ----------------------------------------------------------------------- */
	
	/* Initialized data */
	
	static double c = .564189583547756;
	static double a[5] = { 7.7105849500132e-5,-.00133733772997339,
		.0323076579225834,.0479137145607681,.128379167095513 };
	static double b[3] = { .00301048631703895,.0538971687740286,
		.375795757275549 };
	static double p[8] = { -1.36864857382717e-7,.564195517478974,
		7.21175825088309,43.1622272220567,152.98928504694,
		339.320816734344,451.918953711873,300.459261020162 };
	static double q[8] = { 1.,12.7827273196294,77.0001529352295,
		277.585444743988,638.980264465631,931.35409485061,
		790.950925327898,300.459260956983 };
	static double r[5] = { 2.10144126479064,26.2370141675169,
		21.3688200555087,4.6580782871847,.282094791773523 };
	static double s[4] = { 94.153775055546,187.11481179959,
		99.0191814623914,18.0124575948747 };
	
	/* System generated locals */
	double ret_val, d1;
	
	/* Local variables */
	double e, t, w, ax, bot, top;
	
	/*                     ABS(X) <= 0.5 */
	
	ax = fabs(x);
	if (ax > 0.5) {
		goto L10;
	}
	t = x * x;
	top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.0;
	bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.0;
	ret_val = 0.5 - x * (top / bot) + 0.5;
	if (ind != 0) {
		ret_val = exp(t) * ret_val;
	}
	return ret_val;
	
	/*                  0.5 < ABS(X) <= 4 */
	
L10:
	if (ax > 4.0) {
		goto L20;
	}
	top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax + p[5]) * ax + p[6]) * ax + p[7];
	bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax + q[5]) * ax + q[6]) * ax + q[7];
	ret_val = top / bot;
	goto L40;
	
	/*                      ABS(X) > 4 */
	
L20:
	if (x <= -5.6) {
		goto L50;
	}
	if (ind != 0) {
		goto L30;
	}
	if (x > 100.0) {
		goto L60;
	}
	if (x * x > -exparg(1)) {
		goto L60;
	}
	
L30:
	/* Computing 2nd power */
	d1 = 1.0 / x;
	t = d1 * d1;
	top = (((r[0] * t + r[1]) * t + r[2]) * t + r[3]) * t + r[4];
	bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.0;
	ret_val = (c - t * top / bot) / ax;
	
	/*                      FINAL ASSEMBLY */
	
L40:
	if (ind == 0) {
		goto L41;
	}
	if (x < 0.0) {
		ret_val = exp(x * x) * 2.0 - ret_val;
	}
	return ret_val;
L41:
	w = x * x;
	t = w;
	e = w - t;
	ret_val = (0.5 - e + 0.5) * exp(-t) * ret_val;
	if (x < 0.0) {
		ret_val = 2.0 - ret_val;
	}
	return ret_val;
	
	/*             LIMIT VALUE FOR LARGE NEGATIVE X */
	
L50:
	ret_val = 2.0;
	if (ind != 0) {
		ret_val = exp(x * x) * 2.0;
	}
	return ret_val;
	
	/*             LIMIT VALUE FOR LARGE POSITIVE X */
	/*                       WHEN IND = 0 */
	
L60:
	ret_val = 0.0;
	return ret_val;
}

static double bcorr(double a0, double b0) {
	/* ----------------------------------------------------------------------- */
	/*     EVALUATION OF  DEL(A0) + DEL(B0) - DEL(A0 + B0)  WHERE              */
	/*     LN(GAMMA(A)) = (A - 0.5)*LN(A) - A + 0.5*LN(2*PI) + DEL(A).         */
	/*     IT IS ASSUMED THAT A0 >= 8 AND B0 >= 8.                             */
	/* ----------------------------------------------------------------------- */
	/* Initialized data */
	
	static double c0 = .0833333333333333;
	static double c1 = -.00277777777760991;
	static double c2 = 7.9365066682539e-4;
	static double c3 = -5.9520293135187e-4;
	static double c4 = 8.37308034031215e-4;
	static double c5 = -.00165322962780713;
	
	/* System generated locals */
	double ret_val, r1;
	
	/* Local variables */
	double a, b, c, h, t, w, x, s3, s5, x2, s7, s9, s11;
	/* ----------------------------------------------------------------------- */
	a = RMIN(a0, b0);
	b = RMAX(a0, b0);
	
	h = a / b;
	c = h / (h + 1.0);
	x = 1.0 / (h + 1.0);
	x2 = x * x;
	
	/*                SET SN = (1 - X^N)/(1 - X)                               */
	
	s3 = x + x2 + 1.0;
	s5 = x + x2 * s3 + 1.0;
	s7 = x + x2 * s5 + 1.0;
	s9 = x + x2 * s7 + 1.0;
	s11 = x + x2 * s9 + 1.0;
	
	/*                SET W = DEL(B) - DEL(A + B)                              */
	
	/* Computing 2nd power */
	r1 = 1.0 / b;
	t = r1 * r1;
	w = ((((c5 * s11 * t + c4 * s9) * t + c3 * s7) * t + c2 * s5) * t + c1 * s3) * t + c0;
	w *= c / b;
	
	/*                   COMPUTE  DEL(A) + W                                   */
	
	/* Computing 2nd power */
	r1 = 1.0 / a;
	t = r1 * r1;
	ret_val = (((((c5 * t + c4) * t + c3) * t + c2) * t + c1) * t + c0) / a + w;
	return ret_val;
}

static double basym(double a, double b, double lambda, double eps, int log_p) {
	/* ----------------------------------------------------------------------- */
	/*     ASYMPTOTIC EXPANSION FOR I_x(A,B) FOR LARGE A AND B. */
	/*     LAMBDA = (A + B)*Y - B  AND EPS IS THE TOLERANCE USED. */
	/*     IT IS ASSUMED THAT LAMBDA IS NONNEGATIVE AND THAT */
	/*     A AND B ARE GREATER THAN OR EQUAL TO 15. */
	/* ----------------------------------------------------------------------- */
	/*     ****** NUM IS THE RMAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP */
	/*            ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN. */
	/* ----------------------------------------------------------------------- */
#define num_IT 20
	/*            THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1. */
	
	static double const e0 = 1.12837916709551;/* e0 == 2/sqrt(pi) */
	static double const e1 = .353553390593274;/* e1 == 2^(-3/2)   */
	static double const ln_e0 = 0.120782237635245; /* == ln(e0) */
	
	double a0[num_IT + 1], b0[num_IT + 1], c[num_IT + 1], d[num_IT + 1];
	double f, h, r, s, t, u, w, z, j0, j1, h2, r0, r1, t0, t1, w0, z0, z2, hn, zn;
	double sum, znm1, bsum, dsum;
	
	int i, j, m, n, im1, mm1, np1, imj, mmj;
	
	/* ------------------------ */
	
	f = a * rlog1(-lambda/a) + b * rlog1(lambda/b);
	if(log_p)
		t = -f;
	else {
		t = exp(-f);
		if (t == 0.0) {
			return 0; /* once underflow, always underflow .. */
		}
	}
	z0 = sqrt(f);
	z = z0 / e1 * 0.5;
	z2 = f + f;
	
	if (a < b) {
		h = a / b;
		r0 = 1.0 / (h + 1.0);
		r1 = (b - a) / b;
		w0 = 1.0 / sqrt(a * (h + 1.0));
	} else {
		h = b / a;
		r0 = 1.0 / (h + 1.0);
		r1 = (b - a) / a;
		w0 = 1.0 / sqrt(b * (h + 1.0));
	}
	
	a0[0] = r1 * .66666666666666663;
	c[0] = a0[0] * -0.5;
	d[0] = -c[0];
	j0 = 0.5 / e0 * erfc1(1, z0);
	j1 = e1;
	sum = j0 + d[0] * w0 * j1;
	
	s = 1.0;
	h2 = h * h;
	hn = 1.0;
	w = w0;
	znm1 = z;
	zn = z2;
	for (n = 2; n <= num_IT; n += 2) {
		hn = h2 * hn;
		a0[n - 1] = r0 * 2.0 * (h * hn + 1.0) / (n + 2.0);
		np1 = n + 1;
		s += hn;
		a0[np1 - 1] = r1 * 2.0 * s / (n + 3.0);
		
		for (i = n; i <= np1; ++i) {
			r = (i + 1.0) * -0.5;
			b0[0] = r * a0[0];
			for (m = 2; m <= i; ++m) {
				bsum = 0.0;
				mm1 = m - 1;
				for (j = 1; j <= mm1; ++j) {
					mmj = m - j;
					bsum += (j * r - mmj) * a0[j - 1] * b0[mmj - 1];
				}
				b0[m - 1] = r * a0[m - 1] + bsum / m;
			}
			c[i - 1] = b0[i - 1] / (i + 1.0);
			
			dsum = 0.0;
			im1 = i - 1;
			for (j = 1; j <= im1; ++j) {
				imj = i - j;
				dsum += d[imj - 1] * c[j - 1];
			}
			d[i - 1] = -(dsum + c[i - 1]);
		}
		
		j0 = e1 * znm1 + (n - 1.0) * j0;
		j1 = e1 * zn + n * j1;
		znm1 = z2 * znm1;
		zn = z2 * zn;
		w = w0 * w;
		t0 = d[n - 1] * w * j0;
		w = w0 * w;
		t1 = d[np1 - 1] * w * j1;
		sum += t0 + t1;
		if (fabs(t0) + fabs(t1) <= eps * sum) {
			break;
		}
	}
	
	if(log_p)
		return ln_e0 + t - bcorr(a, b) + log(sum);
	else {
		u = exp(-bcorr(a, b));
		return e0 * t * u * sum;
	}
}

static double alnrel(double a) {
	/* -----------------------------------------------------------------------
	 *            Evaluation of the function ln(1 + a)
	 * ----------------------------------------------------------------------- */
	
	static double p1 = -1.29418923021993;
	static double p2 = .405303492862024;
	static double p3 = -.0178874546012214;
	static double q1 = -1.62752256355323;
	static double q2 = .747811014037616;
	static double q3 = -.0845104217945565;
	
	if (fabs(a) <= 0.375) {
		double t, t2, w;
		t = a / (a + 2.0);
		t2 = t * t;
		w = (((p3 * t2 + p2) * t2 + p1) * t2 + 1.) / (((q3 * t2 + q2) * t2 + q1) * t2 + 1.);
		return t * 2.0 * w;
	} else {
		double x = a + 1.;
		return log(x);
	}
}

static double gamln1(double a) {
	/* ----------------------------------------------------------------------- */
	/*     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 <= A <= 1.25 */
	/* ----------------------------------------------------------------------- */
	
	/* Initialized data */
	
	static double p0 = .577215664901533;
	static double p1 = .844203922187225;
	static double p2 = -.168860593646662;
	static double p3 = -.780427615533591;
	static double p4 = -.402055799310489;
	static double p5 = -.0673562214325671;
	static double p6 = -.00271935708322958;
	static double q1 = 2.88743195473681;
	static double q2 = 3.12755088914843;
	static double q3 = 1.56875193295039;
	static double q4 = .361951990101499;
	static double q5 = .0325038868253937;
	static double q6 = 6.67465618796164e-4;
	static double r0 = .422784335098467;
	static double r1 = .848044614534529;
	static double r2 = .565221050691933;
	static double r3 = .156513060486551;
	static double r4 = .017050248402265;
	static double r5 = 4.97958207639485e-4;
	static double s1 = 1.24313399877507;
	static double s2 = .548042109832463;
	static double s3 = .10155218743983;
	static double s4 = .00713309612391;
	static double s5 = 1.16165475989616e-4;
	
	double w;
	
	if (a < 0.6) {
		w = ((((((p6 * a + p5)* a + p4)* a + p3)* a + p2)* a + p1)* a + p0) / ((((((q6 * a + q5)* a + q4)* a + q3)* a + q2)* a + q1)* a + 1.);
		return -(a) * w;
	}
	else {
		double x = a - 0.5 - 0.5;
		w = (((((r5 * x + r4) * x + r3) * x + r2) * x + r1) * x + r0) / (((((s5 * x + s4) * x + s3) * x + s2) * x + s1) * x + 1.0);
		return x * w;
	}
}

static double gamln(double a) {
	/* ----------------------------------------------------------------------- */
	/*            Evaluation of  ln(gamma(a))  for positive a                  */
	/* ----------------------------------------------------------------------- */
	/*			Written by Alfred H. Morris                                    */
	/*          Naval Surface Warfare Center                                   */
	/*          Dahlgren, Virginia                                             */
	/* ----------------------------------------------------------------------- */
	
	static double d = .418938533204673;/* d == 0.5*(LN(2*PI) - 1) */
	
	static double c0 = .0833333333333333;
	static double c1 = -.00277777777760991;
	static double c2 = 7.9365066682539e-4;
	static double c3 = -5.9520293135187e-4;
	static double c4 = 8.37308034031215e-4;
	static double c5 = -.00165322962780713;
	
	if (a <= 0.8)
		return gamln1(a) - log(a);
	else if (a <= 2.25)
		return gamln1(a - 0.5 - 0.5);
	else if (a < 10.0) {
		int i, n = (int)(a - 1.25);
		double t = a;
		double w = 1.0;
		for (i = 1; i <= n; ++i) {
			t += -1.0;
			w *= t;
		}
		return gamln1(t - 1.) + log(w);
	}
	else { /* a >= 10 */
		double t = 1. / (a * a);
		double w = (((((c5 * t + c4) * t + c3) * t + c2) * t + c1) * t + c0) / a;
		return d + w + (a - 0.5) * (log(a) - 1.0);
	}
}

static double algdiv(double a, double b) {
	/* ----------------------------------------------------------------------- */
	/*     COMPUTATION OF LN(GAMMA(B)/GAMMA(A+B)) WHEN B >= 8                  */
	/*     IN THIS ALGORITHM, DEL(X) IS THE FUNCTION DEFINED BY                */
	/*     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X).         */
	/* ----------------------------------------------------------------------- */
	
	/* Initialized data */
	
	static double c0 = .0833333333333333;
	static double c1 = -.00277777777760991;
	static double c2 = 7.9365066682539e-4;
	static double c3 = -5.9520293135187e-4;
	static double c4 = 8.37308034031215e-4;
	static double c5 = -.00165322962780713;
	
	double c, d, h, t, u, v, w, x, s3, s5, x2, s7, s9, s11;
	
	/* ------------------------ */
	if (a > b) {
		h = b / a;
		c = 1.0 / (h + 1.0);
		x = h / (h + 1.0);
		d = a + (b - 0.5);
	}
	else {
		h = a / b;
		c = h / (h + 1.0);
		x = 1.0 / (h + 1.0);
		d = b + (a - 0.5);
	}
	
	/* Set s<n> = (1 - x^n)/(1 - x) : */
	
	x2 = x * x;
	s3 = x + x2 + 1.0;
	s5 = x + x2 * s3 + 1.0;
	s7 = x + x2 * s5 + 1.0;
	s9 = x + x2 * s7 + 1.0;
	s11 = x + x2 * s9 + 1.0;
	
	/* w := Del(b) - Del(a + b) */
	
	t = 1./ (b * b);
	w = ((((c5 * s11 * t + c4 * s9) * t + c3 * s7) * t + c2 * s5) * t + c1 * s3) * t + c0;
	w *= c / b;
	
	/*                    COMBINE THE RESULTS */
	
	u = d * alnrel(a / b);
	v = a * (log(b) - 1.0);
	if (u > v)
		return w - v - u;
	else
		return w - u - v;
}

static double gsumln(double a, double b) {
	/* ----------------------------------------------------------------------- */
	/*          EVALUATION OF THE FUNCTION LN(GAMMA(A + B))                    */
	/*          FOR 1 <= A <= 2  AND  1 <= B <= 2                              */
	/* ----------------------------------------------------------------------- */
	double x = a + b - 2.;
	
	if (x <= 0.25)
		return gamln1(x + 1.0);
	
	/* else */
	if (x <= 1.25)
		return gamln1(x) + alnrel(x);
	/* else x > 1.25 : */
	return gamln1(x - 1.0) + log(x * (x + 1.0));
}

static double betaln(double a0, double b0) {
	/* -----------------------------------------------------------------------
	 *     Evaluation of the logarithm of the beta function  ln(beta(a0,b0))
	 * ----------------------------------------------------------------------- */
	static double e = .918938533204673;/* e == 0.5*LN(2*PI) */
	
	double a, b, c, h, u, v, w, z;
	int i, n;
	
	a = RMIN(a0 ,b0);
	b = RMAX(a0, b0);
	if (a >= 8.0) {
		goto L60;
	}
	if (a < 1.0) {
		/* ----------------------------------------------------------------------- */
		/*                   PROCEDURE WHEN A < 1                                  */
		/* ----------------------------------------------------------------------- */
		if (b < 8.0)
			return gamln(a) + (gamln(b) - gamln(a+b));
		else
			return gamln(a) + algdiv(a, b);
	}
	/* else */
	/* ----------------------------------------------------------------------- */
	/*                PROCEDURE WHEN 1 <= A < 8                                */
	/* ----------------------------------------------------------------------- */
	if (a > 2.0) {
		goto L30;
	}
	if (b <= 2.0) {
		return gamln(a) + gamln(b) - gsumln(a, b);
	}
	/* else */
	
	w = 0.0;
	if (b < 8.0) {
		goto L40;
	}
	return gamln(a) + algdiv(a, b);
	
L30:
	/*                REDUCTION OF A WHEN B <= 1000                            */
	
	if (b > 1e3) {
		goto L50;
	}
	n = (int)(a - 1.0);
	w = 1.0;
	for (i = 1; i <= n; ++i) {
		a += -1.0;
		h = a / b;
		w *= h / (h + 1.0);
	}
	w = log(w);
	if (b < 8.0) {
		goto L40;
	}
	return w + gamln(a) + algdiv(a, b);
	
L40:
	/*                 REDUCTION OF B WHEN B < 8                               */
	
	n = (int)(b - 1.0);
	z = 1.0;
	for (i = 1; i <= n; ++i) {
		b += -1.0;
		z *= b / (a + b);
	}
	return w + log(z) + (gamln(a) + (gamln(b) - gsumln(a, b)));
	
L50:
	/*                REDUCTION OF A WHEN B > 1000                             */
	n = (int)(a - 1.0);
	w = 1.0;
	for (i = 1; i <= n; ++i) {
		a += -1.0;
		w *= a / (a / b + 1.0);
	}
	return log(w) - n * log(b) + (gamln(a) + algdiv(a, b));
	
L60:
	/* ----------------------------------------------------------------------- */
	/*                   PROCEDURE WHEN A >= 8                                 */
	/* ----------------------------------------------------------------------- */
	
	w = bcorr(a, b);
	h = a / b;
	c = h / (h + 1.0);
	u = -(a - 0.5) * log(c);
	v = b * alnrel(h);
	if (u > v)
		return log(b) * -0.5 + e + w - v - u;
	else
		return log(b) * -0.5 + e + w - u - v;
	
}

static double esum(int mu, double x) {
	/* ----------------------------------------------------------------------- */
	/*                    EVALUATION OF EXP(MU + X) */
	/* ----------------------------------------------------------------------- */
	double w;
	
	if (x > 0.0) {
		goto L10;
	}
	
	if (mu < 0) {
		goto L20;
	}
	w = mu + x;
	if (w > 0.0) {
		goto L20;
	}
	return exp(w);
	
L10:
	if (mu > 0) {
		goto L20;
	}
	w = mu + x;
	if (w < 0.0) {
		goto L20;
	}
	return exp(w);
	
L20:
	w = (double) (mu);
	return exp(w) * exp(x);
}

static double gam1(double a) {
	/*     ------------------------------------------------------------------ */
	/*     COMPUTATION OF 1/GAMMA(A+1) - 1  FOR -0.5 <= A <= 1.5 */
	/*     ------------------------------------------------------------------ */
	
	/* Initialized data */
	
	static double p[7] = { .577215664901533,-.409078193005776,
		-.230975380857675,.0597275330452234,.0076696818164949,
		-.00514889771323592,5.89597428611429e-4 };
	static double q[5] = { 1.,.427569613095214,.158451672430138,
		.0261132021441447,.00423244297896961 };
	static double r[9] = { -.422784335098468,-.771330383816272,
		-.244757765222226,.118378989872749,9.30357293360349e-4,
		-.0118290993445146,.00223047661158249,2.66505979058923e-4,
		-1.32674909766242e-4 };
	static double s1 = .273076135303957;
	static double s2 = .0559398236957378;
	
	double ret_val;
	double d, t, w, bot, top;
	
	t = a;
	d = a - 0.5;
	if (d > 0.0) {
		t = d - 0.5;
	}
	if (t < 0.0) {
		goto L30;
	} else if (t == 0) {
		goto L10;
	} else {
		goto L20;
	}
	
L10:
	ret_val = 0.0;
	return ret_val;
L20:
	top = (((((p[6] * t + p[5]) * t + p[4]) * t + p[3]) * t + p[2]) * t + p[1]) * t + p[0];
	bot = (((q[4] * t + q[3]) * t + q[2]) * t + q[1]) * t + 1.0;
	w = top / bot;
	if (d > 0.0) {
		goto L21;
	}
	ret_val = a * w;
	return ret_val;
L21:
	ret_val = t / a * (w - 0.5 - 0.5);
	return ret_val;
L30:
	top = (((((((r[8] * t + r[7]) * t + r[6]) * t + r[5]) * t + r[4]) * t + r[3]) * t + r[2]) * t + r[1]) * t + r[0];
	bot = (s2 * t + s1) * t + 1.0;
	w = top / bot;
	if (d > 0.0) {
		goto L31;
	}
	ret_val = a * (w + 0.5 + 0.5);
	return ret_val;
L31:
	ret_val = t * w / a;
	return ret_val;
}

static double brcmp1(int mu, double a, double b, double x, double y) {
	/* -----------------------------------------------------------------------
	 *          EVALUATION OF  EXP(MU) * (X^A * Y^B / BETA(A,B))
	 * ----------------------------------------------------------------------- */
	
	static double const__ = .398942280401433; /* == 1/sqrt(2*pi); */
	/* R has  M_1_SQRT_2PI */
	
	/* System generated locals */
	double ret_val, r1;
	
	/* Local variables */
	double c, e, h;
	int i, n;
	double t, u, v, z, a0, b0, x0, y0, apb, lnx, lny;
	double lambda;
	
	a0 = RMIN(a,b);
	if (a0 >= 8.0) {
		goto L100;
	}
	
	if (x > .375) {
		goto L10;
	}
	lnx = log(x);
	lny = alnrel(-x);
	goto L20;
L10:
	if (y > .375) {
		goto L11;
	}
	lnx = alnrel(-y);
	lny = log(y);
	goto L20;
L11:
	lnx = log(x);
	lny = log(y);
L20:
	z = a * lnx + b * lny;
	if (a0 < 1.0) {
		goto L30;
	}
	z -= betaln(a, b);
	ret_val = esum(mu, z);
	return ret_val;
	/* ----------------------------------------------------------------------- */
	/*              PROCEDURE FOR A < 1 OR B < 1 */
	/* ----------------------------------------------------------------------- */
L30:
	b0 = RMAX(a,b);
	if (b0 >= 8.0) {
		goto L80;
	}
	if (b0 > 1.0) {
		goto L60;
	}
	
	/*                   ALGORITHM FOR b0 <= 1 */
	
	ret_val = esum(mu, z);
	if (ret_val == 0.0) {
		return ret_val;
	}
	
	apb = a + b;
	if (apb > 1.0) {
		goto L40;
	}
	z = gam1(apb) + 1.0;
	goto L50;
L40:
	u = a + b - 1.;
	z = (gam1(u) + 1.0) / apb;
	
L50:
	c = (gam1(a) + 1.0) * (gam1(b) + 1.0) / z;
	ret_val = ret_val * (a0 * c) / (a0 / b0 + 1.0);
	return ret_val;
	
	/*                ALGORITHM FOR 1 < b0 < 8 */
	
L60:
	u = gamln1(a0);
	n = (int)(b0 - 1.0);
	if (n < 1) {
		goto L70;
	}
	c = 1.0;
	for (i = 1; i <= n; ++i) {
		b0 += -1.0;
		c *= b0 / (a0 + b0);
		/* L61: */
	}
	u = log(c) + u;
	
L70:
	z -= u;
	b0 += -1.0;
	apb = a0 + b0;
	if (apb > 1.0) {
		goto L71;
	}
	t = gam1(apb) + 1.0;
	goto L72;
L71:
	u = a0 + b0 - 1.;
	t = (gam1(u) + 1.0) / apb;
L72:
	ret_val = a0 * esum(mu, z) * (gam1(b0) + 1.0) / t;
	return ret_val;
	
	/*                   ALGORITHM FOR b0 >= 8 */
	
L80:
	u = gamln1(a0) + algdiv(a0, b0);
	ret_val = a0 * esum(mu, z - u);
	return ret_val;
	/* ----------------------------------------------------------------------- */
	/*              PROCEDURE FOR A >= 8 AND B >= 8 */
	/* ----------------------------------------------------------------------- */
L100:
	if (a > b) {
		goto L101;
	}
	h = a / b;
	x0 = h / (h + 1.0);
	y0 = 1.0 / (h + 1.0);
	lambda = a - (a + b) * x;
	goto L110;
L101:
	h = b / a;
	x0 = 1.0 / (h + 1.0);
	y0 = h / (h + 1.0);
	lambda = (a + b) * y - b;
L110:
	e = -lambda / a;
	if (fabs(e) > 0.6) {
		goto L111;
	}
	u = rlog1(e);
	goto L120;
L111:
	u = e - log(x / x0);
L120:
	e = lambda / b;
	if (fabs(e) > 0.6) {
		goto L121;
	}
	v = rlog1(e);
	goto L130;
L121:
	v = e - log(y / y0);
L130:
	r1 = -(a * u + b * v);
	z = esum(mu, r1);
	
	return const__ * sqrt(b * x0) * z * exp(-bcorr(a, b));
	
}

static double bup(double a, double b, double x, double y, int n, double eps) {
	/* ----------------------------------------------------------------------- */
	/*     EVALUATION OF I_x(A,B) - I_x(A+N,B) WHERE N IS A POSITIVE INT. */
	/*     EPS IS THE TOLERANCE USED. */
	/* ----------------------------------------------------------------------- */
	
	/* System generated locals */
	double ret_val;
	
	/* Local variables */
	int i, k, mu, nm1, kp1;
	double d, l, r, t, w;
	double ap1, apb;
	
	/*          OBTAIN THE SCALING FACTOR EXP(-MU) AND */
	/*             EXP(MU)*(X**A*Y**B/BETA(A,B))/A */
	
	apb = a + b;
	ap1 = a + 1.0;
	mu = 0;
	d = 1.0;
	if (n == 1 || a < 1.0) {
		goto L10;
	}
	if (apb < ap1 * 1.1) {
		goto L10;
	}
	mu = (int)(fabs(exparg(1)));
	k  = (int) exparg(0);
	if (k < mu) {
		mu = k;
	}
	t = (double) mu;
	d = exp(-t);
	
L10:
	ret_val = brcmp1(mu, a, b, x, y) / a;
	if (n == 1 || ret_val == 0.0) {
		return ret_val;
	}
	nm1 = n - 1;
	w = d;
	
	/*          LET K BE THE INDEX OF THE RMAXIMUM TERM */
	
	k = 0;
	if (b <= 1.0) {
		goto L40;
	}
	if (y > 1e-4) {
		goto L20;
	}
	k = nm1;
	goto L30;
L20:
	r = (b - 1.0) * x / y - a;
	if (r < 1.0) {
		goto L40;
	}
	k = nm1;
	t = (double) nm1;
	if (r < t) {
		k = (int) r;
	}
	
	/*          ADD THE INCREASING TERMS OF THE SERIES */
	
L30:
	for (i = 1; i <= k; ++i) {
		l = (double) (i - 1);
		d = (apb + l) / (ap1 + l) * x * d;
		w += d;
		/* L31: */
	}
	if (k == nm1) {
		goto L50;
	}
	
	/*          ADD THE REMAINING TERMS OF THE SERIES */
	
L40:
	kp1 = k + 1;
	for (i = kp1; i <= nm1; ++i) {
		l = (double) (i - 1);
		d = (apb + l) / (ap1 + l) * x * d;
		w += d;
		if (d <= eps * w) {
			goto L50;
		}
	}
	
	/*               TERRMINATE THE PROCEDURE */
	
L50:
	ret_val *= w;
	return ret_val;
}

static double psi(double x) {
	/* ---------------------------------------------------------------------
	 
	 *                 Evaluation of the Digamma function psi(x)
	 
	 *                           -----------
	 
	 *     Psi(xx) is assigned the value 0 when the digamma function cannot
	 *     be computed.
	 
	 *     The main computation involves evaluation of rational Chebyshev
	 *     approximations published in Math. Comp. 27, 123-127(1973) by
	 *     Cody, Strecok and Thacher.                                        */
	/* --------------------------------------------------------------------- */
	/*     Psi was written at Argonne National Laboratory for the FUNPACK */
	/*     package of special function subroutines. Psi was modified by */
	/*     A.H. Morris (NSWC). */
	/* --------------------------------------------------------------------- */
	static double piov4 = .785398163397448; /* == pi / 4 */
	/*     dx0 = zero of psi() to extended precision :                       */
	static double dx0 = 1.461632144968362341262659542325721325;
	/* --------------------------------------------------------------------- */
	/*     COEFFICIENTS FOR RATIONAL APPROXIMATION OF                        */
	/*     PSI(X) / (X - X0),  0.5 <= X <= 3.0                               */
	static double p1[7] = { .0089538502298197,4.77762828042627,
		142.441585084029,1186.45200713425,3633.51846806499,
		4138.10161269013,1305.60269827897 };
	static double q1[6] = { 44.8452573429826,520.752771467162,
		2210.0079924783,3641.27349079381,1908.310765963,
		6.91091682714533e-6 };
	/* --------------------------------------------------------------------- */
	/*     COEFFICIENTS FOR RATIONAL APPROXIMATION OF                        */
	/*     PSI(X) - LN(X) + 1 / (2*X),  X > 3.0                              */
	static double p2[4] = { -2.12940445131011,-7.01677227766759,
		-4.48616543918019,-.648157123766197 };
	static double q2[4] = { 32.2703493791143,89.2920700481861,
		54.6117738103215,7.77788548522962 };
	/* --------------------------------------------------------------------- */
	int i, m, n, nq;
	double d2;
	double w, z;
	double den, aug, sgn, xmx0, xRMAX1, upper, xsmall;
	/* --------------------------------------------------------------------- */
	/*     MACHINE DEPENDENT CONSTANTS ...                                   */
	/* --------------------------------------------------------------------- */
	/*	  XRMAX1	 = THE SMALLEST POSITIVE FLOATING POINT CONSTANT
	 WITH ENTIRELY INT REPRESENTATION.  ALSO USED
	 AS NEGATIVE OF LOWER BOUND ON ACCEPTABLE NEGATIVE
	 ARGUMENTS AND AS THE POSITIVE ARGUMENT BEYOND WHICH
	 PSI MAY BE REPRESENTED AS LOG(X).
	 * Originally:  xRMAX1 = aRMIN1(ipmpar(3), 1./spmpar(1))                 */
	//xRMAX1 = (double) INT_RMAX;
	xRMAX1 = Rf_d1mach(4) - 1.0; // Daniel: changed on 03/01/2009 according to R-2.4.0, R-2.6.0
	d2 = 0.5 / Rf_d1mach(3);
	if(xRMAX1 > d2) xRMAX1 = d2;
	
	/* --------------------------------------------------------------------- */
	/*        XSMALL = ABSOLUTE ARGUMENT BELOW WHICH PI*COTAN(PI*X)          */
	/*                 MAY BE REPRESENTED BY 1/X.                            */
	xsmall = 1e-9;
	/* --------------------------------------------------------------------- */
	aug = 0.0;
	if (x < 0.5) {
		/* --------------------------------------------------------------------- */
		/*     X < 0.5,  USE REFLECTION FORMULA                                  */
		/*     PSI(1-X) = PSI(X) + PI * COTAN(PI*X)                              */
		/* --------------------------------------------------------------------- */
		if (fabs(x) <= xsmall) {
			if (x == 0.0) {
				goto L_err;
			}
			/* --------------------------------------------------------------------- */
			/*     0 < ABS(X) <= XSMALL.  USE 1/X AS A SUBSTITUTE                    */
			/*     FOR  PI*COTAN(PI*X)                                               */
			/* --------------------------------------------------------------------- */
			aug = -1.0 / x;
		} else { /* |x| > xsmall */
			/* --------------------------------------------------------------------- */
			/*     REDUCTION OF ARGUMENT FOR COTAN                                   */
			/* --------------------------------------------------------------------- */
			/* L100: */
			w = -x;
			sgn = piov4;
			if (w <= 0.0) {
				w = -w;
				sgn = -sgn;
			}
			/* --------------------------------------------------------------------- */
			/*     MAKE AN ERROR EXIT IF |X| >= XRMAX1                               */
			/* --------------------------------------------------------------------- */
			if (w >= xRMAX1) {
				goto L_err;
			}
			nq = (int) w;
			w -= (double) nq;
			nq = (int) (w * 4.0);
			w = (w - (double) nq * 0.25) * 4.0;
			/* --------------------------------------------------------------------- */
			/*     W IS NOW RELATED TO THE FRACTIONAL PART OF  4.0 * X.              */
			/*     ADJUST ARGUMENT TO CORRESPOND TO VALUES IN FIRST                  */
			/*     QUADRANT AND DETERRMINE SIGN                                      */
			/* --------------------------------------------------------------------- */
			n = nq / 2;
			if (n + n != nq) {
				w = 1.0 - w;
			}
			z = piov4 * w;
			m = n / 2;
			if (m + m != n) {
				sgn = -sgn;
			}
			/* --------------------------------------------------------------------- */
			/*     DETERRMINE FINAL VALUE FOR  -PI*COTAN(PI*X)                       */
			/* --------------------------------------------------------------------- */
			n = (nq + 1) / 2;
			m = n / 2;
			m += m;
			if (m == n) {
				/* --------------------------------------------------------------------- */
				/*     CHECK FOR SINGULARITY                                             */
				/* --------------------------------------------------------------------- */
				if (z == 0.0) {
					goto L_err;
				}
				/* --------------------------------------------------------------------- */
				/*     USE COS/SIN AS A SUBSTITUTE FOR COTAN, AND                        */
				/*     SIN/COS AS A SUBSTITUTE FOR TAN                                   */
				/* --------------------------------------------------------------------- */
				aug = sgn * (cos(z) / sin(z) * 4.0);
			} else { /* L140: */
				aug = sgn * (sin(z) / cos(z) * 4.0);
			}
		}
		x = 1.0 - x;
	}
	/* L200: */
	if (x <= 3.0) {
		/* --------------------------------------------------------------------- */
		/*     0.5 <= X <= 3.0                                                   */
		/* --------------------------------------------------------------------- */
		den = x;
		upper = p1[0] * x;
		
		for (i = 1; i <= 5; ++i) {
			den = (den + q1[i - 1]) * x;
			upper = (upper + p1[i]) * x;
		}
		
		den = (upper + p1[6]) / (den + q1[5]);
		xmx0 = x - dx0;
		return den * xmx0 + aug;
	}
	/* --------------------------------------------------------------------- */
	/*     IF X >= XRMAX1, PSI = LN(X)                                       */
	/* --------------------------------------------------------------------- */
	if (x < xRMAX1) {
		/* --------------------------------------------------------------------- */
		/*     3.0 < X < XRMAX1                                                  */
		/* --------------------------------------------------------------------- */
		w = 1.0 / (x * x);
		den = w;
		upper = p2[0] * w;
		
		for (i = 1; i <= 3; ++i) {
			den = (den + q2[i - 1]) * w;
			upper = (upper + p2[i]) * w;
		}
		
		aug = upper / (den + q2[3]) - 0.5 / x + aug;
	}
	return aug + log(x);
	/* --------------------------------------------------------------------- */
	/*     ERROR RETURN                                                      */
	/* --------------------------------------------------------------------- */
L_err:
	return 0.;
}

double fpser(double a, double b, double x, double eps, int log_p) {
	/* ----------------------------------------------------------------------- *
	 
	 *                 EVALUATION OF I (A,B)
	 *                                X
	 
	 *          FOR B < RMIN(EPS, EPS*A) AND X <= 0.5
	 
	 * ----------------------------------------------------------------------- */
	
	double ans, c, s, t, an, tol;
	
	/* SET  ans := x^a : */
	if (log_p) {
		ans = a * log(x);
	}
	else if (a > eps * 0.001) {
		t = a * log(x);
		if (t < exparg(1)) { /* exp(t) would underflow */
			return 0.0;
		}
		ans = exp(t);
	}
	else
		ans = 1.;
	
	/*                NOTE THAT 1/B(A,B) = B */
	
	if (log_p)
		ans += log(b) - log(a);
	else
		ans *= b / a;
	
	tol = eps / a;
	an = a + 1.0;
	t = x;
	s = t / an;
	do {
		an += 1.0;
		t   = x * t;
		c   = t / an;
		s  += c;
	} while (fabs(c) > tol);
	
	if (log_p)
		ans += r_log1p(a * s);
	else
		ans *= a * s + 1.0;
	return ans;
}

static double apser(double a, double b, double x, double eps) {
	/* -----------------------------------------------------------------------
	 *     apser() yields the incomplete beta ratio  I_{1-x}(b,a)  for
	 *     a <= RMIN(eps,eps*b), b*x <= 1, and x <= 0.5,  i.e., a is very small.
	 *     Use only if above inequalities are satisfied.
	 * ----------------------------------------------------------------------- */
	
	static double const g = .577215664901533;
	
	double tol, c, j, s, t, aj;
	double bx = b * x;
	
	t = x - bx;
	if (b * eps <= 0.02)
		c = log(x) + psi(b) + g + t;
	else
		c = log(bx) + g + t;
	
	tol = eps * 5.0 * fabs(c);
	j = 1.;
	s = 0.;
	do {
		j += 1.0;
		t *= x - bx / j;
		aj = t / j;
		s += aj;
	} while (fabs(aj) > tol);
	
	return -a * (c + s);
}

static double bpser(double a, double b, double x, double eps, int log_p) {
	/* -----------------------------------------------------------------------
	 * Power SERies expansion for evaluating I_x(a,b) when
	 *	       b <= 1 or b*x <= 0.7.   eps is the tolerance used.
	 * ----------------------------------------------------------------------- */
	
	int i, m;
	double ans, c, n, t, u, w, z, a0, b0, apb, tol, sum;
	
	if (x == 0.) {
		return R_D__0;
	}
	/* ----------------------------------------------------------------------- */
	/*	      compute the factor  x^a/(a*Beta(a,b)) */
	/* ----------------------------------------------------------------------- */
	a0 = RMIN(a,b);
	if (a0 >= 1.0) { /*		 ------	 1 <= a0 <= b0  ------ */
		z = a * log(x) - betaln(a, b);
		ans = log_p ? z - log(a) : exp(z) / a;
	}
	else {
		b0 = RMAX(a,b);
		
		if (b0 < 8.0) {
			if (b0 <= 1.0) { /*	 ------	 a0 < 1	 and  b0 <= 1  ------ */
				if(log_p) {
					ans = a * log(x);
				}
				else {
					ans = pow(x, a);
					if (ans == 0.) /* once underflow, always underflow .. */
						return ans;
				}
				apb = a + b;
				if (apb > 1.0) {
					u = a + b - 1.;
					z = (gam1(u) + 1.0) / apb;
				}
				else {
					z = gam1(apb) + 1.0;
				}
				c = (gam1(a) + 1.0) * (gam1(b) + 1.0) / z;
				
				if(log_p) /* FIXME ? -- improve quite a bit for c ~= 1 */
					ans += log(c * (b / apb));
				else
					ans *=  c * (b / apb);
			}
			else { /* 	------	a0 < 1 < b0 < 8	 ------ */
				u = gamln1(a0);
				m = (int)(b0 - 1.0);
				if (m >= 1) {
					c = 1.0;
					for (i = 1; i <= m; ++i) {
						b0 += -1.0;
						c  *= b0 / (a0 + b0);
					}
					u += log(c);
				}
				
				z = a * log(x) - u;
				b0 += -1.0;
				apb = a0 + b0;
				if (apb > 1.0) {
					u = a0 + b0 - 1.;
					t = (gam1(u) + 1.0) / apb;
				}
				else {
					t = gam1(apb) + 1.0;
				}
				
				if(log_p) /* FIXME? potential for improving log(t) */
					ans = z + log(a0 / a) + r_log1p(gam1(b0)) - log(t);
				else
					ans = exp(z) * (a0 / a) * (gam1(b0) + 1.0) / t;
			}
			
		}
		else { /* 		------  a0 < 1 < 8 <= b0  ------ */
			u = gamln1(a0) + algdiv(a0, b0);
			z = a * log(x) - u;
			
			if(log_p)
				ans = z + log(a0 / a);
			else
				ans = a0 / a * exp(z);
		}
	}
	
	if (!log_p && (ans == 0.0 || a <= eps * 0.1)) {
		return ans;
	}
	
	/* ----------------------------------------------------------------------- */
	/*		       COMPUTE THE SERIES */
	/* ----------------------------------------------------------------------- */
	sum = 0.;
	n = 0.;
	c = 1.;
	tol = eps / a;
	
	do {
		n += 1.;
		c *= (0.5 - b / n + 0.5) * x;
		w = c / (a + n);
		sum += w;
	} while (fabs(w) > tol);
	
	if(log_p)
		ans += r_log1p(a * sum);
	else
		ans *= a * sum + 1.0;
	return ans;
}

static double brcomp(double a, double b, double x, double y, int log_p) {
	/* -----------------------------------------------------------------------
	 *		 Evaluation of x^a * y^b / Beta(a,b)
	 * ----------------------------------------------------------------------- */
	
	static double const__ = .398942280401433; /* == 1/sqrt(2*pi); */
	/* R has  M_1_SQRT_2PI , and M_LN_SQRT_2PI = ln(sqrt(2*pi)) = 0.918938.. */
	
	int i, n;
	double c, e, h, t, u, v, z, a0, b0, x0, y0, apb, lnx, lny;
	double lambda;
	
	
	if (x == 0.0 || y == 0.0) {
		return R_D__0;
	}
	a0 = RMIN(a, b);
	if (a0 >= 8.0) {
		goto L100;
	}
	
	if (x <= .375) {
		lnx = log(x);
		lny = alnrel(-x);
	}
	else {
		if (y > .375) {
			lnx = log(x);
			lny = log(y);
		} else {
			lnx = alnrel(-y);
			lny = log(y);
		}
	}
	
	z = a * lnx + b * lny;
	if (a0 >= 1.) {
		z -= betaln(a, b);
		return R_D_exp(z);
	}
	
	/* ----------------------------------------------------------------------- */
	/*		PROCEDURE FOR a < 1 OR b < 1 */
	/* ----------------------------------------------------------------------- */
	
	b0 = RMAX(a, b);
	if (b0 >= 8.0) { /* L80: */
		u = gamln1(a0) + algdiv(a0, b0);
		return (log_p ? log(a0) + (z - u) : a0 * exp(z - u));
	}
	/* else : */
	
	if (b0 <= 1.0) { /*		algorithm for RMAX(a,b) = b0 <= 1 */
		double e_z = R_D_exp(z);
		
		if (!log_p && e_z == 0.0) /* exp() underflow */
			return 0.;
		
		apb = a + b;
		if (apb > 1.0) {
			u = a + b - 1.;
			z = (gam1(u) + 1.0) / apb;
		} else {
			z = gam1(apb) + 1.0;
		}
		
		c = (gam1(a) + 1.0) * (gam1(b) + 1.0) / z;
		/* FIXME? log(a0*c)= log(a0)+ log(c) and that is improvable */
		return (log_p ? e_z + log(a0 * c) - r_log1p(a0/b0) : e_z * (a0 * c) / (a0 / b0 + 1.0));
	}
	/* else : */
	
	/*		  ALGORITHM FOR 1 < b0 < 8 */
	
	u = gamln1(a0);
	n = (int)(b0 - 1.0);
	if (n >= 1) {
		c = 1.0;
		for (i = 1; i <= n; ++i) {
			b0 += -1.0;
			c *= b0 / (a0 + b0);
		}
		u = log(c) + u;
	}
	z -= u;
	b0 += -1.0;
	apb = a0 + b0;
	if (apb > 1.0) {
		u = a0 + b0 - 1.;
		t = (gam1(u) + 1.0) / apb;
	} else {
		t = gam1(apb) + 1.0;
	}
	
	return (log_p ? log(a0) + z + r_log1p(gam1(b0)) - log(t) : a0 * exp(z) * (gam1(b0) + 1.0) / t);
	
	/* ----------------------------------------------------------------------- */
	/*		PROCEDURE FOR A >= 8 AND B >= 8 */
	/* ----------------------------------------------------------------------- */
L100:
	if (a <= b) {
		h  = a / b;
		x0 = h / (h + 1.0);
		y0 = 1.0 / (h + 1.0);
		lambda = a - (a + b) * x;
	} else {
		h = b / a;
		x0 = 1.0 / (h + 1.0);
		y0 = h / (h + 1.0);
		lambda = (a + b) * y - b;
	}
	
	e = -lambda / a;
	if (fabs(e) > .6)
		u = e - log(x / x0);
	else
		u = rlog1(e);
	
	e = lambda / b;
	if (fabs(e) <= .6)
		v = rlog1(e);
	else
		v = e - log(y / y0);
	
	z = log_p ? -(a * u + b * v) : exp(-(a * u + b * v));
	
	return(log_p ? -M_LN_SQRT_2PI + .5*log(b * x0) + z - bcorr(a,b) : const__ * sqrt(b * x0) * z * exp(-bcorr(a, b)));
}

static double bfrac(double a, double b, double x, double y, double lambda, double eps, int log_p) {
	/* -----------------------------------------------------------------------
	 Continued fraction expansion for I_x(a,b) when a, b > 1.
	 It is assumed that  lambda = (a + b)*y - b.
	 -----------------------------------------------------------------------*/
	
	double c, e, n, p, r, s, t, w, c0, c1, r0, an, bn, yp1, anp1, bnp1, beta, alpha;
	
	double brc = brcomp(a, b, x, y, log_p);
	
	if (!log_p && brc == 0.) /* already underflowed to 0 */
		return 0.;
	
	c   = lambda + 1.0;
	c0  = b / a;
	c1  = 1.0 / a + 1.0;
	yp1 = y + 1.0;
	
	n    = 0.0;
	p    = 1.0;
	s    = a + 1.0;
	an   = 0.0;
	bn   = 1.0;
	anp1 = 1.0;
	bnp1 = c / c1;
	r    = c1 / c;
	
	/*        CONTINUED FRACTION CALCULATION */
	
	do {
		n += 1.0;
		t = n / a;
		w = n * (b - n) * x;
		e = a / s;
		alpha = p * (p + c0) * e * e * (w * x);
		e = (t + 1.0) / (c1 + t + t);
		beta = n + w / s + e * (c + n * yp1);
		p = t + 1.0;
		s += 2.0;
		
		/* update an, bn, anp1, and bnp1 */
		
		t = alpha * an + beta * anp1;
		an = anp1;
		anp1 = t;
		t = alpha * bn + beta * bnp1;
		bn = bnp1;
		bnp1 = t;
		
		r0 = r;
		r = anp1 / bnp1;
		if (fabs(r - r0) <= eps * r) {
			break;
		}
		
		/* rescale an, bn, anp1, and bnp1 */
		
		an /= bnp1;
		bn /= bnp1;
		anp1 = r;
		bnp1 = 1.0;
	} while (1);
	
	return (log_p ? brc + log(r) : brc * r);
}

static double erf__(double x) {
	/* -----------------------------------------------------------------------
	 *             EVALUATION OF THE REAL ERROR FUNCTION
	 * ----------------------------------------------------------------------- */
	
	/* Initialized data */
	
	static double c = .564189583547756;
	static double a[5] = { 7.7105849500132e-5,-.00133733772997339,
		.0323076579225834,.0479137145607681,.128379167095513 };
	static double b[3] = { .00301048631703895,.0538971687740286,
		.375795757275549 };
	static double p[8] = { -1.36864857382717e-7,.564195517478974,
		7.21175825088309,43.1622272220567,152.98928504694,
		339.320816734344,451.918953711873,300.459261020162 };
	static double q[8] = { 1.,12.7827273196294,77.0001529352295,
		277.585444743988,638.980264465631,931.35409485061,
		790.950925327898,300.459260956983 };
	static double r[5] = { 2.10144126479064,26.2370141675169,
		21.3688200555087,4.6580782871847,.282094791773523 };
	static double s[4] = { 94.153775055546,187.11481179959,
		99.0191814623914,18.0124575948747 };
	
	/* System generated locals */
	double ret_val;
	
	/* Local variables */
	double t, x2, ax, bot, top;
	
	ax = fabs(x);
	if (ax <= 0.5) {
		t = x * x;
		top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.0;
		bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.0;
		
		return x * (top / bot);
	}
	/* else: ax > 0.5 */
	
	if (ax <= 4.) { /*  ax in (0.5, 4] */
		top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax + p[5]) * ax + p[6]) * ax + p[7];
		bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax + q[5]) * ax + q[6]) * ax + q[7];
		ret_val = 0.5 - exp(-x * x) * top / bot + 0.5;
		if (x < 0.0) {
			ret_val = -ret_val;
		}
		return ret_val;
	}
	
	/* else: ax > 4 */
	
	if (ax >= 5.8) {
		return x > 0 ? 1 : -1;
	}
	x2 = x * x;
	t = 1.0 / x2;
	top = (((r[0] * t + r[1]) * t + r[2]) * t + r[3]) * t + r[4];
	bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.0;
	t = (c - top / (x2 * bot)) / ax;
	ret_val = 0.5 - exp(-x2) * t + 0.5;
	if (x < 0.0) {
		ret_val = -ret_val;
	}
	return ret_val;
}

static void grat1(double a, double x, double r, double *p, double *q, double eps) {
	/* -----------------------------------------------------------------------
	 *        Evaluation of the incomplete gamma ratio functions
	 *                      P(a,x) and Q(a,x)
	 
	 *     It is assumed that a <= 1.  eps is the tolerance to be used.
	 *     the input argument r has the value  r = e^(-x)* x^a / Gamma(a).
	 * ----------------------------------------------------------------------- */
	
	double c, g, h, j, l, t, w, z, an, am0, an0, a2n, b2n, cma;
	double tol, sum, a2nm1, b2nm1;
	
	if (a * x == 0.0) { /* L130: */
		if (x <= a)
			goto L100;
		else
			goto L110;
	}
	else if (a == 0.5) {
		goto L120;
	}
	
	if (x < 1.1) { /* L10:  Taylor series for  P(a,x)/x^a */
		an = 3.0;
		c = x;
		sum = x / (a + 3.0);
		tol = eps * 0.1 / (a + 1.0);
		do {
			an += 1.0;
			c = -c * (x / an);
			t = c / (a + an);
			sum += t;
		} while (fabs(t) > tol);
		
		j = a * x * ((sum / 6.0 - 0.5 / (a + 2.0)) * x + 1.0 / (a + 1.0));
		
		z = a * log(x);
		h = gam1(a);
		g = h + 1.0;
		if (x >= 0.25) {
			if (a < x / 2.59) {
				goto L40;
			}
		}
		else {
			if (z > -0.13394) {
				goto L40;
			}
		}
		
		w = exp(z);
		*p = w * g * (0.5 - j + 0.5);
		*q = 0.5 - *p + 0.5;
		return;
		
	L40:
		l = rr_expm1(z);
		w = l + 0.5 + 0.5;
		*q = (w * j - l) * g - h;
		if (*q < 0.0) {
			goto L110;
		}
		*p = 0.5 - *q + 0.5;
		return;
	}
	
	/* L50: ----  (x >= 1.1)  ---- Continued Fraction Expansion */
	
	a2nm1 = 1.0;
	a2n = 1.0;
	b2nm1 = x;
	b2n = x + (1.0 - a);
	c = 1.0;
	
	do {
		a2nm1 = x * a2n + c * a2nm1;
		b2nm1 = x * b2n + c * b2nm1;
		am0 = a2nm1 / b2nm1;
		c += 1.0;
		cma = c - a;
		a2n = a2nm1 + cma * a2n;
		b2n = b2nm1 + cma * b2n;
		an0 = a2n / b2n;
	} while (fabs(an0 - am0) >= eps * an0);
	
	*q = r * an0;
	*p = 0.5 - *q + 0.5;
	return;
	
	/*                SPECIAL CASES */
	
L100:
	*p = 0.0;
	*q = 1.0;
	return;
L110:
	*p = 1.0;
	*q = 0.0;
	return;
L120:
	if (x < 0.25) {
		*p = erf__(sqrt(x));
		*q = 0.5 - *p + 0.5;
	} else {
		*q = erfc1(0, sqrt(x));
		*p = 0.5 - *q + 0.5;
	}
	return;
	
}

static void bgrat(double a, double b, double x, double y, double *w, double eps, int *ierr) {
	/* -----------------------------------------------------------------------
	 *     Asymptotic Expansion for I_x(A,B)  when a is larger than b.
	 *     The result of the expansion is added to w.
	 *     It is assumed a >= 15 and b <= 1.
	 *     eps is the tolerance used.
	 *     ierr is a variable that reports the status of the results.
	 * ----------------------------------------------------------------------- */
	
	double c[30], d[30];
	int i, n, nm1;
	double j, l, p, q, r, s, t, u, v, z, n2, t2, dj, cn, nu, bm1;
	double lnx, sum, bp2n, coef;
	
	bm1 = b - 0.5 - 0.5;
	nu = a + bm1 * 0.5;
	if (y > 0.375)
		lnx = log(x);
	else
		lnx = alnrel(-y);
	
	z = -nu * lnx;
	if (b * z == 0.0) {
		goto L_Error;
	}
	
	/*                 COMPUTATION OF THE EXPANSION */
	
	/* set r := exp(-z) * z^b / Gamma(b) */
	r = b * (gam1(b) + 1.0) * exp(b * log(z));
	
	r = r * exp(a * lnx) * exp(bm1 * 0.5 * lnx);
	u = algdiv(b, a) + b * log(nu);
	u = r * exp(-u);
	if (u == 0.0) {
		goto L_Error;
	}
	grat1(b, z, r, &p, &q, eps); /* -> (p,q)  {p + q = 1} */
	
	v = 0.25 / (nu * nu);
	t2 = lnx * 0.25 * lnx;
	l = *w / u;
	j = q / r;
	sum = j;
	t = 1.0;
	cn = 1.0;
	n2 = 0.0;
	for (n = 1; n <= 30; ++n) {
		bp2n = b + n2;
		j = (bp2n * (bp2n + 1.0) * j + (z + bp2n + 1.0) * t) * v;
		n2 += 2.0;
		t *= t2;
		cn /= n2 * (n2 + 1.0);
		nm1 = n - 1;
		c[nm1] = cn;
		s = 0.0;
		if (n > 1) {
			coef = b - n;
			for (i = 1; i <= nm1; ++i) {
				s += coef * c[i - 1] * d[nm1 - i];
				coef += b;
			}
		}
		d[nm1] = bm1 * cn + s / n;
		dj = d[nm1] * j;
		sum += dj;
		if (sum <= 0.0) {
			goto L_Error;
		}
		if (fabs(dj) <= eps * (sum + l)) {
			break;
		}
	}
	
	/*                    ADD THE RESULTS TO W */
	*ierr = 0;
	*w += u * sum;
	return;
	
	/*               THE EXPANSION CANNOT BE COMPUTED */
	
L_Error:
	*ierr = 1;
	return;
}

void bratio(double a, double b, double x, double y, double *w, double *w1, int *ierr, int log_p) {
	/* -----------------------------------------------------------------------
	 
	 *	      Evaluation of the Incomplete Beta function I_x(a,b)
	 
	 *		       --------------------
	 
	 *     It is assumed that a and b are nonnegative, and that x <= 1
	 *     and y = 1 - x.  Bratio assigns w and w1 the values
	 
	 *			w  = I_x(a,b)
	 *			w1 = 1 - I_x(a,b)
	 
	 *     ierr is a variable that reports the status of the results.
	 *     If no input errors are detected then ierr is set to 0 and
	 *     w and w1 are computed. otherwise, if an error is detected,
	 *     then w and w1 are assigned the value 0 and ierr is set to
	 *     one of the following values ...
	 
	 *	  ierr = 1  if a or b is negative
	 *	  ierr = 2  if a = b = 0
	 *	  ierr = 3  if x < 0 or x > 1
	 *	  ierr = 4  if y < 0 or y > 1
	 *	  ierr = 5  if x + y != 1
	 *	  ierr = 6  if x = a = 0
	 *	  ierr = 7  if y = b = 0
	 
	 * --------------------
	 *     Written by Alfred H. Morris, Jr.
	 *	  Naval Surface Warfare Center
	 *	  Dahlgren, Virginia
	 *     Revised ... Nov 1991
	 * ----------------------------------------------------------------------- */
	
	bool do_swap;
	int n, ierr1;
	double z, a0, b0, x0, y0, eps, lambda;
	
	/*  eps is a machine dependent constant: the smallest floating point number for which   1.0 + eps > 1.0 */
	eps = 2.0 * Rf_d1mach(3); /* == DBL_EPSILON (in R, Rmath) */
	
	/* ----------------------------------------------------------------------- */
	*w  = R_D__0;
	*w1 = R_D__0;
	
	if (a < 0.0 || b < 0.0)   { *ierr = 1; return; }
	if (a == 0.0 && b == 0.0) { *ierr = 2; return; }
	if (x < 0.0 || x > 1.0)   {	*ierr = 3; return; }
	if (y < 0.0 || y > 1.0)   { *ierr = 4; return; }
	
	z = x + y - 0.5 - 0.5;
	
	if (fabs(z) > eps * 3.0) { *ierr = 5; return; }
	
	*ierr = 0;
	if (x == 0.0) goto L200;
	if (y == 0.0) goto L210;
	
	if (a == 0.0) goto L211;
	if (b == 0.0) goto L201;
	
	eps = RMAX(eps, 1e-15);
	if (RMAX(a,b) < eps * .001) { /* procedure for a and b < 0.001 * eps */
		/* L230: */
		if(log_p) {
			z   = log(a + b);
			*w	= log(b) - z;
			*w1 = log(a) - z;
		} else {
			*w	= b / (a + b);
			*w1 = a / (a + b);
		}
		return;
	}
	
#define SET_0_noswap \
	a0 = a;  x0 = x; \
	b0 = b;  y0 = y;
	
#define SET_0_swap   \
	a0 = b;  x0 = y; \
	b0 = a;  y0 = x;
	
	if (RMIN(a,b) > 1.0) {
		goto L30;
	}
	
	/*             PROCEDURE FOR a0 <= 1 OR b0 <= 1 */
	
	do_swap = (x > 0.5);
	if (do_swap) {
		SET_0_swap;
	} else {
		SET_0_noswap;
	}
	/* now have  x0 <= 1/2 <= y0  (still  x0+y0 == 1) */
	
	if (b0 < RMIN(eps, eps * a0)) { /* L80: */
		*w  = fpser(a0, b0, x0, eps, log_p);
		*w1 = log_p ? R_D_LExp_toms708(*w) : 0.5 - *w + 0.5;
		goto L_end_after_log;
	}
	
	if (a0 < RMIN(eps, eps * b0) && b0 * x0 <= 1.0) { /* L90: */
		*w1 = apser(a0, b0, x0, eps);
		*w  = 0.5 - *w1 + 0.5;
		goto L_end;
	}
	
	if (RMAX(a0,b0) > 1.0) { /* L20:  RMIN(a,b) <= 1 < RMAX(a,b)  */
		if (b0 <= 1.0) {
			goto L100;
		}
		if (x0 >= 0.3) {
			goto L110;
		}
		if (x0 < 0.1) {
			if (pow(x0*b0, a0) <= 0.7) {
				goto L100;
			}
		}
		if (b0 > 15.0) {
			*w1 = 0.;
			goto L131;
		}
	} else { /*  a, b <= 1 */
		if (a0 >= RMIN(0.2, b0)) {
			goto L100;
		}
		if (pow(x0, a0) <= 0.9) {
			goto L100;
		}
		if (x0 >= 0.3) {
			goto L110;
		}
	}
	n = 20;
	goto L130;
	
L30:
	/*             PROCEDURE FOR a0 > 1 AND b0 > 1 */
	if (a > b)
		lambda = (a + b) * y - b;
	else
		lambda = a - (a + b) * x;
	
	do_swap = (lambda < 0.0);
	if (do_swap) {
		lambda = -lambda;
		SET_0_swap;
	} else {
		SET_0_noswap;
	}
	
	if (b0 < 40.0) {
		if (b0 * x0 <= 0.7)
			goto L100;
		else
			goto L140;
	}
	else if (a0 > b0) { /* ----  a0 > b0 >= 40  ---- */
		if (b0 <= 100.0) {
			goto L120;
		}
		if (lambda > b0 * 0.03) {
			goto L120;
		}
	}
	else if (a0 <= 100.0) {
		goto L120;
	}
	else if (lambda > a0 * 0.03) {
		goto L120;
	}
	
	/* else if none of the above    L180: */
	*w  = basym(a0, b0, lambda, eps * 100.0, log_p);
	*w1 = log_p ? R_D_LExp_toms708(*w) : 0.5 - *w + 0.5;
	goto L_end_after_log;
	
	/*            EVALUATION OF THE APPROPRIATE ALGORITHM */
	
L100:
	*w  = bpser(a0, b0, x0, eps, log_p);
	*w1 = log_p ? R_D_LExp_toms708(*w) : 0.5 - *w + 0.5;
	goto L_end_after_log;
	
L110:
	*w1 = bpser(b0, a0, y0, eps, log_p);
	*w  = log_p ? R_D_LExp_toms708(*w1) : 0.5 - *w1 + 0.5;
	goto L_end_after_log;
	
L120:
	*w = bfrac(a0, b0, x0, y0, lambda, eps * 15.0, log_p);
	*w1 = log_p ? R_D_LExp_toms708(*w) : 0.5 - *w + 0.5;
	goto L_end_after_log;
	
L130:
	*w1 = bup(b0, a0, y0, x0, n, eps);
	b0 += n;
L131:
	bgrat(b0, a0, y0, x0, w1, 15*eps, &ierr1);
	*w = 0.5 - *w1 + 0.5;
	goto L_end;
	
L140:
	/* b0 := fractional_part( b0 )  in (0, 1]  */
	n   = (int) b0;
	b0 -= n;
	if (b0 == 0.) {
		--n;
		b0 = 1.;
	}
	
	*w = bup(b0, a0, y0, x0, n, eps);
	if (x0 <= 0.7) {
		/* log_p :  TODO:  w = bup(.) + bpser(.)  -- not so easy to use log-scale */
		*w += bpser(a0, b0, x0, eps, /* log_p = */ FALSE);
		*w1 = 0.5 - *w + 0.5;
		goto L_end;
	}
	/* L150: */
	if (a0 <= 15.0) {
		n   = 20;
		*w += bup(a0, b0, x0, y0, n, eps);
		a0 += n;
	}
	bgrat(a0, b0, x0, y0, w, 15*eps, &ierr1);
	*w1 = 0.5 - *w + 0.5;
	goto L_end;
	
	
	/* TERRMINATION OF THE PROCEDURE */
	
L200:
	if (a == 0.0) {
		*ierr = 6;
		return;
	}
L201:
	*w  = R_D__0;
	*w1 = R_D__1;
	return;
L210:
	if (b == 0.0) {
		*ierr = 7;
		return;
	}
L211:
	*w  = R_D__1;
	*w1 = R_D__0;
	return;
L_end:
	if (log_p) {
		*w  = log(*w);
		*w1 = log(*w1);
	}
L_end_after_log:
	if (do_swap) { /* swap */
		double t = *w;
		*w = *w1;
		*w1 = t;
	}
	return;
#undef SET_0_noswap
#undef SET_0_swap
}

/*  DESCRIPTION
 *		Compute the log gamma correction factor for x >= 10 so that log(gamma(x)) = .5*log(2*pi) + (x-.5)*log(x) -x + lgammacor(x)
 *		[ lgammacor(x) is called Del(x)	in other contexts (e.g. dcdflib)]
 *  NOTES
 *		This routine is a translation into C of a Fortran subroutine written by W. Fullerton of Los Alamos Scientific Laboratory.
 *  SEE ALSO
 *		Loader(1999)'s stirlerr() {in ./stirlerr.c} is *very* similar in spirit, is faster and cleaner, but is only defined "fast" for half integers.
 */
double lgammacor(double x) {
	const static double algmcs[15] = {
		+.1666389480451863247205729650822e+0,
		-.1384948176067563840732986059135e-4,
		+.9810825646924729426157171547487e-8,
		-.1809129475572494194263306266719e-10,
		+.6221098041892605227126015543416e-13,
		-.3399615005417721944303330599666e-15,
		+.2683181998482698748957538846666e-17,
		-.2868042435334643284144622399999e-19,
		+.3962837061046434803679306666666e-21,
		-.6831888753985766870111999999999e-23,
		+.1429227355942498147573333333333e-24,
		-.3547598158101070547199999999999e-26,
		+.1025680058010470912000000000000e-27,
		-.3401102254316748799999999999999e-29,
		+.1276642195630062933333333333333e-30
	};
	
	double tmp;
	
	/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
	 *   xbig = 2 ^ 26.5
	 *   xmax = DBL_MAX / 48 =  2^1020 / 3 */
#define nalgm 5
#define xbig_lgammacor  94906265.62425156
#define xmax_lgammacor  3.745194030963158e306
	
	if (x < 10)
		ML_ERR_return_NAN
		else if (x >= xmax_lgammacor) {
			ML_ERROR(ME_UNDERFLOW, "lgammacor"); /* allow to underflow below */
		}
		else if (x < xbig_lgammacor) {
			tmp = 10 / x;
			return chebyshev_eval(tmp * tmp * 2 - 1, algmcs, nalgm) / x;
		}
	return 1 / (x * 12);
}

/*  DESCRIPTION
 *		The function lgammafn computes log|gamma(x)|.  The function lgammafn_sign in addition assigns the sign of the gamma function to the address in the second argument if this is not NULL.
 *  NOTES
 *		This routine is a translation into C of a Fortran subroutine by W. Fullerton of Los Alamos Scientific Laboratory.
 *		The accuracy of this routine compares (very) favourably with those of the Sun Microsystems portable mathematical library.
 */
double lgammafn_sign(double x, int *sgn) {
	double ans, y, sinpiy;
	
#ifdef NOMORE_FOR_THREADS
	static double xmax = 0.;
	static double dxrel = 0.;
	
	if (xmax == 0) {/* initialize machine dependent constants _ONCE_ */
		xmax = d1mach(2)/log(d1mach(2));/* = 2.533 e305	 for IEEE double */
		dxrel = sqrt (d1mach(4));/* sqrt(Eps) ~ 1.49 e-8  for IEEE double */
	}
#else
	/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
	 xmax  = DBL_MAX / log(DBL_MAX) = 2^1024 / (1024 * log(2)) = 2^1014 / log(2)
	 dxrel = sqrt(DBL_EPSILON) = 2^-26 = 5^26 * 1e-26 (is *exact* below !)
	 */
#define xmax_lgammafn_sign  2.5327372760800758e+305
#define dxrel               1.490116119384765696e-8
#endif
	
	if (sgn != NULL)
		*sgn = 1;
	
	if(ISNAN(x)!=0)
		return x;
	
	if (x < 0 && fmod(floor(-x), 2.) == 0) //Daniel: I am not sure if there are something missed here.
		
		if (sgn != NULL)
			*sgn = -1;
	
	if (x <= 0 && x == ftrunc(x)) { /* Negative integer argument */
		ML_ERROR(ME_RANGE, "lgamma");
		return ML_POSINF;/* +Inf, since lgamma(x) = log|gamma(x)| */
	}
	
	y = fabs(x);
	
	if (y <= 10)
		return log(fabs(gammafn(x)));
	
	/* ELSE  y = |x| > 10 ---------------------- */
	
	if (y > xmax_lgammafn_sign) {
		ML_ERROR(ME_RANGE, "lgamma");
		return ML_POSINF;
	}
	
	if (x > 0) { /* i.e. y = x > 10 */
		if(x > 1e17)
			return(x*(log(x) - 1.));
		else if(x > 4934720.)
			return(M_LN_SQRT_2PI + (x - 0.5) * log(x) - x);
		else
			return M_LN_SQRT_2PI + (x - 0.5) * log(x) - x + lgammacor(x);
	}
	/* else: x < -10; y = -x */
	sinpiy = fabs(sin(MY_PI * y));
	
	if (sinpiy == 0) { /* Negative integer argument === Now UNNECESSARY: caught above */
		MATHLIB_WARNING(" ** should NEVER happen! *** [lgamma.c: Neg.int, y=%g]\n",y);
		ML_ERR_return_NAN;
	}
	
	ans = M_LN_SQRT_PId2 + (x - 0.5) * log(y) - x - log(sinpiy) - lgammacor(y);
	
	if(fabs((x - ftrunc(x - 0.5)) * ans / x) < dxrel) {
		
		/* The answer is less than half precision because the argument is too near a negative integer. */
		
		ML_ERROR(ME_PRECISION, "lgamma");
	}
	
	return ans;
}
double lgammafn(double x) {
	return lgammafn_sign(x, NULL);
}

/*  DESCRIPTION
 *		Computes the log of the error term in Stirling's formula.
 *      For n > 15, uses the series 1/12n - 1/360n^3 + ...
 *      For n <=15, integers or half-integers, uses stored values.
 *      For other n < 15, uses lgamma directly (don't use this to write lgamma!)
 *
 *  Merge in to R:
 *  Copyright (C) 2000, The R Core Development Team
 *  R has lgammafn, and lgamma is not part of ISO C
 *
 *  stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n )
 *              = log Gamma(n+1) - 1/2 * [log(2*pi) + log(n)] - n*[log(n) - 1]
 *              = log Gamma(n+1) - (n + 1/2) * log(n) + n - log(2*pi)/2
 *  see also lgammacor() in ./lgammacor.c  which computes almost the same!
 */
double stirlerr(double n) {
	
#define S0 0.083333333333333333333       /* 1/12 */
#define S1 0.00277777777777777777778     /* 1/360 */
#define S2 0.00079365079365079365079365  /* 1/1260 */
#define S3 0.000595238095238095238095238 /* 1/1680 */
#define S4 0.0008417508417508417508417508/* 1/1188 */
	
	/*
  error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
	 */
	const static double sferr_halves[31] = {
		0.0, /* n=0 - wrong, place holder only */
		0.1534264097200273452913848,  /* 0.5 */
		0.0810614667953272582196702,  /* 1.0 */
		0.0548141210519176538961390,  /* 1.5 */
		0.0413406959554092940938221,  /* 2.0 */
		0.03316287351993628748511048, /* 2.5 */
		0.02767792568499833914878929, /* 3.0 */
		0.02374616365629749597132920, /* 3.5 */
		0.02079067210376509311152277, /* 4.0 */
		0.01848845053267318523077934, /* 4.5 */
		0.01664469118982119216319487, /* 5.0 */
		0.01513497322191737887351255, /* 5.5 */
		0.01387612882307074799874573, /* 6.0 */
		0.01281046524292022692424986, /* 6.5 */
		0.01189670994589177009505572, /* 7.0 */
		0.01110455975820691732662991, /* 7.5 */
		0.010411265261972096497478567, /* 8.0 */
		0.009799416126158803298389475, /* 8.5 */
		0.009255462182712732917728637, /* 9.0 */
		0.008768700134139385462952823, /* 9.5 */
		0.008330563433362871256469318, /* 10.0 */
		0.007934114564314020547248100, /* 10.5 */
		0.007573675487951840794972024, /* 11.0 */
		0.007244554301320383179543912, /* 11.5 */
		0.006942840107209529865664152, /* 12.0 */
		0.006665247032707682442354394, /* 12.5 */
		0.006408994188004207068439631, /* 13.0 */
		0.006171712263039457647532867, /* 13.5 */
		0.005951370112758847735624416, /* 14.0 */
		0.005746216513010115682023589, /* 14.5 */
		0.005554733551962801371038690  /* 15.0 */
	};
	double nn;
	
	if (n <= 15.0) {
		nn = n + n;
		if (nn == (int)nn)
			return(sferr_halves[(int)nn]);
		return(lgammafn(n + 1.) - (n + 0.5)*log(n) + n - M_LN_SQRT_2PI);
	}
	
	nn = n*n;
	if (n>500) return((S0-S1/nn)/n);
	if (n> 80) return((S0-(S1-S2/nn)/nn)/n);
	if (n> 35) return((S0-(S1-(S2-S3/nn)/nn)/nn)/n);
	/* 15 < n <= 35 : */
	return((S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n);
}

/*  DESCRIPTION
 *		This function computes the value of the gamma function.
 *  NOTES
 *		This function is a translation into C of a Fortran subroutine by W. Fullerton of Los Alamos Scientific Laboratory. (e.g. http://www.netlib.org/slatec/fnlib/gamma.f)
 *		The accuracy of this routine compares (very) favourably with those of the Sun Microsystems portable mathematical library.
 *		MM specialized the case of  n!  for n < 50 - for even better precision
 */
double gammafn(double x) {
	const static double gamcs[42] = {
		+.8571195590989331421920062399942e-2,
		+.4415381324841006757191315771652e-2,
		+.5685043681599363378632664588789e-1,
		-.4219835396418560501012500186624e-2,
		+.1326808181212460220584006796352e-2,
		-.1893024529798880432523947023886e-3,
		+.3606925327441245256578082217225e-4,
		-.6056761904460864218485548290365e-5,
		+.1055829546302283344731823509093e-5,
		-.1811967365542384048291855891166e-6,
		+.3117724964715322277790254593169e-7,
		-.5354219639019687140874081024347e-8,
		+.9193275519859588946887786825940e-9,
		-.1577941280288339761767423273953e-9,
		+.2707980622934954543266540433089e-10,
		-.4646818653825730144081661058933e-11,
		+.7973350192007419656460767175359e-12,
		-.1368078209830916025799499172309e-12,
		+.2347319486563800657233471771688e-13,
		-.4027432614949066932766570534699e-14,
		+.6910051747372100912138336975257e-15,
		-.1185584500221992907052387126192e-15,
		+.2034148542496373955201026051932e-16,
		-.3490054341717405849274012949108e-17,
		+.5987993856485305567135051066026e-18,
		-.1027378057872228074490069778431e-18,
		+.1762702816060529824942759660748e-19,
		-.3024320653735306260958772112042e-20,
		+.5188914660218397839717833550506e-21,
		-.8902770842456576692449251601066e-22,
		+.1527474068493342602274596891306e-22,
		-.2620731256187362900257328332799e-23,
		+.4496464047830538670331046570666e-24,
		-.7714712731336877911703901525333e-25,
		+.1323635453126044036486572714666e-25,
		-.2270999412942928816702313813333e-26,
		+.3896418998003991449320816639999e-27,
		-.6685198115125953327792127999999e-28,
		+.1146998663140024384347613866666e-28,
		-.1967938586345134677295103999999e-29,
		+.3376448816585338090334890666666e-30,
		-.5793070335782135784625493333333e-31
	};
	
	int i, n;
	double y;
	double sinpiy, value;
	
#ifdef NOMORE_FOR_THREADS
	static int ngam = 0;
	static double xmin = 0, xmax = 0., xsml = 0., dxrel = 0.;
	
	/* Initialize machine dependent constants, the first time gamma() is called.
	 FIXME for threads ! */
	if (ngam == 0) {
		ngam = chebyshev_init(gamcs, 42, DBL_EPSILON/20);/*was .1*d1mach(3)*/
		gammalims(&xmin, &xmax);/*-> ./gammalims.c */
		xsml = exp(fmax2(log(DBL_MIN), -log(DBL_MAX)) + 0.01);
		/*   = exp(.01)*DBL_MIN = 2.247e-308 for IEEE */
		dxrel = sqrt(DBL_EPSILON);/*was sqrt(d1mach(4)) */
	}
#else
	/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
	 * (xmin, xmax) are non-trivial, see ./gammalims.c
	 * xsml = exp(.01)*DBL_MIN
	 * dxrel = sqrt(DBL_EPSILON) = 2 ^ -26
	 */
# define ngam 22
# define xmin -170.5674972726612
# define xmax  171.61447887182298
# define xsml 2.2474362225598545e-308
# define dxrel 1.490116119384765696e-8
#endif
	
	if(ISNAN(x)!=0)
		return x;
	
	/* If the argument is exactly zero or a negative integer then return NaN. */
	if (x == 0 || (x < 0 && x == (long)x)) {
		ML_ERROR(ME_DOMAIN, "gammafn");
		return ML_NAN;
	}
	
	y = fabs(x);
	
	if (y <= 10) {
		/* Compute gamma(x) for -10 <= x <= 10
		 * Reduce the interval and find gamma(1 + y) for 0 <= y < 1
		 * first of all. */
		
		n = (int)x;
		if(x < 0)
			--n;
		y = x - n;/* n = floor(x)  ==>	y in [ 0, 1 ) */
		--n;
		value = chebyshev_eval(y * 2 - 1, gamcs, ngam) + .9375;
		if (n == 0)
			return value;/* x = 1.dddd = 1+y */
		
		if (n < 0) {
			/* compute gamma(x) for -10 <= x < 1 */
			/* exact 0 or "-n" checked already above */
			/* The answer is less than half precision because x too near a negative integer. */
			if (x < -0.5 && fabs(x - (int)(x - 0.5) / x) < dxrel) {
				ML_ERROR(ME_PRECISION, "gammafn");
			}
			
			/* The argument is so close to 0 that the result would overflow. */
			if (y < xsml) {
				ML_ERROR(ME_RANGE, "gammafn");
				if(x > 0)
					return ML_POSINF;
				else
					return ML_NEGINF;
			}
			
			n = -n;
			
			for (i = 0; i < n; i++) {
				value /= (x + i);
			}
			return value;
		}
		else {
			/* gamma(x) for 2 <= x <= 10 */
			
			for (i = 1; i <= n; i++) {
				value *= (y + i);
			}
			return value;
		}
	}
	else {
		/* gamma(x) for	 y = |x| > 10. */
		
		if (x > xmax) {			/* Overflow */
			ML_ERROR(ME_RANGE, "gammafn");
			return ML_POSINF;
		}
		
		if (x < xmin) {			/* Underflow */
			ML_ERROR(ME_UNDERFLOW, "gammafn");
			return 0.;
		}
		
		if(y <= 50 && y == (int)y) { /* compute (n - 1)! */
			value = 1.;
			for (i = 2; i < y; i++)
				value *= i;
		}
		else { /* normal case */
			value = exp((y - 0.5) * log(y) - y + M_LN_SQRT_2PI + ((2*y == (int)2*y)? stirlerr(y) : lgammacor(y)));
		}
		if (x > 0)
			return value;
		
		if (fabs((x - (int)(x - 0.5))/x) < dxrel){
			/* The answer is less than half precision because the argument is too near a negative integer. */
			ML_ERROR(ME_PRECISION, "gammafn");
		}
		
		sinpiy = sin(MY_PI * y);
		if (sinpiy == 0) {		/* Negative integer arg - overflow */
			ML_ERROR(ME_RANGE, "gammafn");
			return ML_POSINF;
		}
		
		return -MY_PI / (y * sinpiy * value);
	}
}

/* DESCRIPTION
 *		To compute the binomial probability, call dbinom(x,n,p). This checks for argument validity, and calls dbinom_raw().
 *		dbinom_raw() does the actual computation; note this is called by other functions in addition to dbinom().
 *			(1) dbinom_raw() has both p and q arguments, when one may be represented more accurately than the other (in particular, in df()).
 *			(2) dbinom_raw() does NOT check that inputs x and n are integers. This should be done in the calling function, where necessary.
 *				-- but is not the case at all when called e.g., from df() or dbeta() !
 *			(3) Also does not check for 0 <= p <= 1 and 0 <= q <= 1 or NaN's. Do this in the calling function.
 */
double dbinom_raw(double x, double n, double p, double q, int give_log) {
	double lf, lc;
	
	if (p == 0)
		return((x == 0) ? R_D__1 : R_D__0);
	if (q == 0)
		return((x == n) ? R_D__1 : R_D__0);
	
	if (x == 0) {
		if(n == 0)
			return R_D__1;
		lc = (p < 0.1) ? -bd0(n,n*q) - n*p : n*log(q);
		return( R_D_exp(lc) );
	}
	if (x == n) {
		lc = (q < 0.1) ? -bd0(n,n*p) - n*q : n*log(p);
		return( R_D_exp(lc) );
	}
	if (x < 0 || x > n)
		return( R_D__0 );
	
	/* n*p or n*q can underflow to zero if n and p or q are small.  This used to occur in dbeta, and gives NaN as from R 2.3.0.  */
	lc = stirlerr(n) - stirlerr(x) - stirlerr(n-x) - bd0(x,n*p) - bd0(n-x,n*q);
	
	/* f = (M_2PI*x*(n-x))/n; could overflow or underflow */
	/* Upto R 2.7.1:
	 * lf = log(M_2PI) + log(x) + log(n-x) - log(n);
	 * -- following is much better for  x << n : */
	lf = log(M_2PI) + log(x) + r_log1p(- x/n);
	
	return R_D_exp(lc - 0.5*lf);
}

/*  DESCRIPTION
 *		This function returns the value of the log beta function.
 *  NOTES
 *		This routine is a translation into C of a Fortran subroutine by W. Fullerton of Los Alamos Scientific Laboratory.
 */
double lbeta(double a, double b) {
	double corr, p, q;
	
	p = q = a;
	if(b < p) p = b;/* := min(a,b) */
	if(b > q) q = b;/* := max(a,b) */
	
	if(ISNAN(a)!=0 || ISNAN(b)!=0)
		return a + b;
	
	/* both arguments must be >= 0 */
	
	if (p < 0)
		ML_ERR_return_NAN
		else if (p == 0) {
			return ML_POSINF;
		}
		else if (FINITE(q)==0) {
			return ML_NEGINF;
		}
	
	if (p >= 10) {
		/* p and q are big. */
		corr = lgammacor(p) + lgammacor(q) - lgammacor(p + q);
		return log(q) * -0.5 + M_LN_SQRT_2PI + corr + (p - 0.5) * log(p / (p + q)) + q * r_log1p(-p / (p + q));
	}
	else if (q >= 10) {
		/* p is small, but q is big. */
		corr = lgammacor(q) - lgammacor(p + q);
		return lgammafn(p) + corr + p - p * log(p + q) + (q - 0.5) * r_log1p(-p / (p + q));
	}
	else
	/* p and q are small: p <= q < 10. */
		return log(gammafn(p) * (gammafn(q) / gammafn(p + q)));
}

/*  DESCRIPTION
 *		Returns distribution function of the beta distribution. ( = The incomplete beta ratio I_x(p,q) ).
 *  NOTES
 *		As from R 2.3.0, a wrapper for TOMS708
 *      as from R 2.6.0, 'log_p' partially improved over log(p..)
 */
double pbeta_raw(double x, double pin, double qin, int lower_tail, int log_p) {
	double x1 = 0.5 - x + 0.5, w, wc;
	int ierr;
	bratio(pin, qin, x, x1, &w, &wc, &ierr, log_p); /* -> ./toms708.c */
	if(ierr)
		MATHLIB_WARNING("pbeta_raw() -> bratio() gave error code %d", ierr);
	return lower_tail ? w : wc;
}

/*  DESCRIPTION
 *		the quantile function of Beta distribution
 *		qbeta(double alpha, double p, double q, int lower_tail, int log_p)
 *			alpha     : probability / quantile
 *			p         : shape parameter
 *			q         : scale parameter
 *			lower_tail: should be set to 1
 *			log_p     : shoule be set to 0
 *  REFERENCE
 *		Cran, G. W., K. J. Martin and G. E. Thomas (1977). Remark AS R19 and Algorithm AS 109, Applied Statistics, 26(1), 111-114.
 *		Remark AS R83 (v.39, 309-310) and the correction (v.40(1) p.236) have been incorporated in this version.
 *
 *	set the exponent of accu to -2r-2 for r digits of accuracy
 *	---- NEW ---- -- still fails for p = 1e11, q=.5
 */
#define fpu     3e-308
#define acu_min 1e-300 /* acu_min:  Minimal value for accuracy 'acu' which will depend on (a,p); acu_min >= fpu ! */
#define lower   fpu
#define upper   1-2.22e-16

#define const1 2.30753
#define const2 0.27061
#define const3 0.99229
#define const4 0.04481

double qbeta(double alpha, double p, double q, int lower_tail, int log_p) {
	int swap_tail, i_pb, i_inn;
	double a, adj, logbeta, g, h, pp, p_, prev, qq, r, s, t, tx, w, y, yprev;
	double acu;
	volatile double xinbta;
	
	/* test for admissibility of parameters */
	
	if (ISNAN(p)!=0 || ISNAN(q)!=0 || ISNAN(alpha)!=0)
		return p + q + alpha;
	if(p < 0. || q < 0.)
		ML_ERR_return_NAN;
	
	R_Q_P01_boundaries(alpha, 0, 1);
	
	p_ = R_DT_qIv(alpha);/* lower_tail prob (in any case) */
	
	if(log_p && (p_ == 0. || p_ == 1.))
		return p_; /* better than NaN or infinite loop; FIXME: suboptimal, since -Inf < alpha ! */
	
	/* initialize */
	logbeta = lbeta(p, q);
	
	/* change tail if necessary;  afterwards   0 < a <= 1/2	 */
	if (p_ <= 0.5) {
		a = p_;
		pp = p;
		qq = q;
		swap_tail = 0;
	}
	else { /* change tail, swap  p <-> q :*/
		a = (!lower_tail && !log_p)? alpha : 1 - p_;
		pp = q;
		qq = p;
		swap_tail = 1;
	}
	
	/* calculate the initial approximation */
	
	/* y := {fast approximation of} qnorm(1 - a) :*/
	r = sqrt(-2 * log(a));
	y = r - (const1 + const2 * r) / (1. + (const3 + const4 * r) * r);
	if (pp > 1 && qq > 1) {
		r = (y * y - 3.) / 6.;
		s = 1. / (pp + pp - 1.);
		t = 1. / (qq + qq - 1.);
		h = 2. / (s + t);
		w = y * sqrt(h + r) / h - (t - s) * (r + 5. / 6. - 2. / (3. * h));
		xinbta = pp / (pp + qq * exp(w + w));
	}
	else {
		r = qq + qq;
		t = 1. / (9. * qq);
		t = r * pow(1. - t + y * sqrt(t), 3.0);
		if (t <= 0.)
			xinbta = 1. - exp((r_log1p(-a)+ log(qq) + logbeta) / qq);
		else {
			t = (4. * pp + r - 2.) / t;
			if (t <= 1.)
				xinbta = exp((log(a * pp) + logbeta) / pp);
			else
				xinbta = 1. - 2. / (t + 1.);
		}
	}
	
	/* solve for x by a modified newton-raphson method, using the function pbeta_raw */
	
	r = 1 - pp;
	t = 1 - qq;
	yprev = 0.;
	adj = 1;
	/* Sometimes the approximation is negative! */
	if (xinbta < lower)
		xinbta = 0.5;
	else if (xinbta > upper)
		xinbta = 0.5;
	
	/* Desired accuracy should depend on  (a,p)
	 * This is from Remark .. on AS 109, adapted.
	 * However, it's not clear if this is "optimal" for IEEE double prec.
	 * acu = fmax2(acu_min, pow(10., -25. - 5./(pp * pp) - 1./(a * a)));
	 * NEW: 'acu' accuracy NOT for squared adjustment, but simple;
	 * ---- i.e.,  "new acu" = sqrt(old acu)
	 */
	acu = fmax2(acu_min, pow(10., -13 - 2.5/(pp * pp) - 0.5/(a * a)));
	tx = prev = 0.;	/* keep -Wall happy */
	
	for (i_pb=0; i_pb < 1000; i_pb++) {
		y = pbeta_raw(xinbta, pp, qq, /*lower_tail = */ TRUE, FALSE);
		
		if(FINITE(y)==0)
			ML_ERR_return_NAN;
		
		y = (y - a) * exp(logbeta + r * log(xinbta) + t * r_log1p(-xinbta));
		if (y * yprev <= 0.)
			prev = fmax2(fabs(adj),fpu);
		g = 1;
		for (i_inn=0; i_inn < 1000;i_inn++) {
			adj = g * y;
			if (fabs(adj) < prev) {
				tx = xinbta - adj; /* trial new x */
				if (tx >= 0. && tx <= 1) {
					if (prev <= acu)	      goto L_converged;
					if (fabs(y) <= acu)       goto L_converged;
					if (tx != 0. && tx != 1)  break;
				}
			}
			g /= 3;
		}
		if (fabs(tx - xinbta) < 1e-15*xinbta) goto L_converged;
		xinbta = tx;
		yprev = y;
	}
	/*-- NOT converged: Iteration count --*/
	ML_ERROR(ME_PRECISION, "qbeta");
	
L_converged:
	return swap_tail ? 1 - xinbta : xinbta;
}

/*  DESCRIPTION
 *    Beta density,
 *                   (a+b-1)!     a-1       b-1
 *      p(x;a,b) = ------------ x     (1-x)
 *                 (a-1)!(b-1)!
 *
 *               = (a+b-1) dbinom(a-1; a+b-2,x)
 *
 *    The basic formula for the log density is thus (a-1) log x + (b-1) log (1-x) - lbeta(a, b)
 *    If either a or b <= 2 then 0 < lbeta(a, b) < 710 and so no term is large.  We use Loader's code only if both a and b > 2.
 */
double dbeta(double x, double a, double b, int give_log) {
	double lval;
	
	if (ISNAN(x)!=0 || ISNAN(a)!=0 || ISNAN(b)!=0)
		return x + a + b;
	
	if (a <= 0 || b <= 0)
		ML_ERR_return_NAN;
	if (x < 0 || x > 1)
		return(R_D__0);
	if (x == 0) {
		if(a > 1)
			return(R_D__0);
		if(a < 1)
			return(ML_POSINF);
		/* a == 1 : */
		return(R_D_val(b));
	}
	if (x == 1) {
		if(b > 1)
			return(R_D__0);
		if(b < 1)
			return(ML_POSINF);
		/* b == 1 : */
		return(R_D_val(a));
	}
	if (a <= 2 || b <= 2)
		lval = (a-1)*log(x) + (b-1)*r_log1p(-x) - lbeta(a, b);
	else
		lval = log(a+b-1) + dbinom_raw(a-1, a+b-2, x, 1-x, TRUE);
	
	return R_D_exp(lval);
}

//================================================================================================= Likelihood and Data
void ResetParams() {
	weight.clear();
	condLike.clear();
	pMatrix.assign(16, 0.0);
	baseFreq.assign(4, 0.25);
	kappa  = 1.0;
	brlens = BRLENS_PR_MEAN;
	lnlike = 0.0;
#if PHYML
	unscaledBaseFreq.assign(4, 25.0);
	brlensMin = 1.E-8;
	brlensMax = 100.0;
	min_diff_lk_local = 1.E-04;
	optimizeBaseFreq  = false;
#endif	
}

void ReformatData(std::string &s, std::vector<int> &q, std::vector<int> &weight, std::vector<double> &condLike) {
	// reformat either tumor or matched normal pileup
	//
	typedef std::map<int, int>::value_type valType;
	std::map<int, int> AQualityCount;
	std::map<int, int> CQualityCount;
	std::map<int, int> GQualityCount;
	std::map<int, int> TQualityCount;
	unsigned missingCount = 0;
	for(unsigned i = 0; i < s.length(); i++) {
		int quality = q[i];
		switch(s[i]) {
			case '?':
				missingCount += 1;
				break;
			case 'A': {
				std::map<int, int>::iterator itA = AQualityCount.find(quality);
				if(itA != AQualityCount.end())
					itA->second += 1;
				else
					AQualityCount.insert(valType(quality, 1));
			}
				break;
			case 'C': {
				std::map<int, int>::iterator itC = CQualityCount.find(quality);
				if(itC != CQualityCount.end())
					itC->second += 1;
				else
					CQualityCount.insert(valType(quality, 1));
			}
				break;
			case 'G': {
				std::map<int, int>::iterator itG = GQualityCount.find(quality);
				if(itG != GQualityCount.end())
					itG->second += 1;
				else
					GQualityCount.insert(valType(quality, 1));
			}
				break;
			case 'T': {
				std::map<int, int>::iterator itT = TQualityCount.find(quality);
				if(itT != TQualityCount.end())
					itT->second += 1;
				else
					TQualityCount.insert(valType(quality, 1));
			}
				break;
			default:
				std::cout << "Unknown character in sequences.\n";
				exit(-1);
				break;
		}
	}
	
	unsigned totalPattern = unsigned(AQualityCount.size() + CQualityCount.size() + GQualityCount.size() + TQualityCount.size());
	if(missingCount != 0)
		totalPattern += 1;
	weight.assign(totalPattern, 0);
	condLike.assign(totalPattern*4, 0.0);
	std::vector<int>::iterator    itWeight   = weight.begin();
	std::vector<double>::iterator itCondLike = condLike.begin();
	
	for(std::map<int, int>::iterator it = AQualityCount.begin(); it != AQualityCount.end(); it++) {
		*itWeight = it->second;
		itWeight++;
		double probIncorrect = pow(10.0, (double)(it->first)/(-10.0));
		*itCondLike       = 1.0 - probIncorrect;
		*(itCondLike + 1) = probIncorrect/3.0;
		*(itCondLike + 2) = probIncorrect/3.0;
		*(itCondLike + 3) = probIncorrect/3.0;
		itCondLike += 4;
	}
	for(std::map<int, int>::iterator it = CQualityCount.begin(); it != CQualityCount.end(); it++) {
		*itWeight = it->second;
		itWeight++;
		double probIncorrect = pow(10.0, (double)(it->first)/(-10.0));
		*itCondLike       = probIncorrect/3.0;
		*(itCondLike + 1) = 1.0 - probIncorrect;
		*(itCondLike + 2) = probIncorrect/3.0;
		*(itCondLike + 3) = probIncorrect/3.0;
		itCondLike += 4;
	}
	for(std::map<int, int>::iterator it = GQualityCount.begin(); it != GQualityCount.end(); it++) {
		*itWeight = it->second;
		itWeight++;
		double probIncorrect = pow(10.0, (double)(it->first)/(-10.0));
		*itCondLike       = probIncorrect/3.0;
		*(itCondLike + 1) = probIncorrect/3.0;
		*(itCondLike + 2) = 1.0 - probIncorrect;
		*(itCondLike + 3) = probIncorrect/3.0;
		itCondLike += 4;
	}
	for(std::map<int, int>::iterator it = TQualityCount.begin(); it != TQualityCount.end(); it++) {
		*itWeight = it->second;
		itWeight++;
		double probIncorrect = pow(10.0, (double)(it->first)/(-10.0));
		*itCondLike       = probIncorrect/3.0;
		*(itCondLike + 1) = probIncorrect/3.0;
		*(itCondLike + 2) = probIncorrect/3.0;
		*(itCondLike + 3) = 1.0 - probIncorrect;
		itCondLike += 4;
	}
	if(missingCount != 0) {
		*itWeight = missingCount;
		*itCondLike       = 1.0;
		*(itCondLike + 1) = 1.0;
		*(itCondLike + 2) = 1.0;
		*(itCondLike + 3) = 1.0;
	}
}

void CalcProbMat(std::vector<double> &pMat, const std::vector<double> &state_freqs, const double &kappa, const double &brlens) {
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

#if PHYML

#define ROUND_MAX      100
#define UNLIKELY       -1.e10
#define BRENT_IT_MAX   500
#define BRENT_ZEPS     1.e-10
#define BRENT_CGOLD    0.3819660
#define SMALL          DBL_MIN
#define YES            1
#define NO             0
#define ALF            1.0e-4
#define ITMAX          200
#define EPS            3.0e-8
#define STPMX          100.0
#define FABS           fabs
#define FLOOR          floor
#define POW            pow
#define SQRT           sqrt
#define SIGN(a,b)      ((b) > 0.0 ? fabs(a) : -fabs(a))
#define For(i,n)       for(i=0; i<n; i++)
#define SHFT(a,b,c,d)  (a)=(b);(b)=(c);(c)=(d);
#define MAX(a,b)       ((a)>(b)?(a):(b))

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) < SMALL ? 0.0 : sqrarg*sqrarg)


void *mCalloc(int nb, size_t size) {
	void *allocated;
	if((allocated = calloc((size_t)nb,size)) != NULL) {
		return allocated;
    }
	else {
		std::cerr << "\nError: low memory\n";
		exit(1);
	}
	return NULL;
}

void *mRealloc(void *p,int nb, size_t size) {
	if((p = realloc(p,(size_t)nb*size)) != NULL)
		return p;
	else {
		std::cerr << "\nError: low memory\n";
		exit(1);
	}
	return NULL;
}

void Free(void *p) {
	free(p);
	p = NULL;
}

void Check_Br_Len_Bounds() {
	if(brlens > brlensMax) brlens = brlensMax;
	if(brlens < brlensMin) brlens = brlensMin;
}

void NormalizeBaseFreq() {
	double sum;
	int i;
	
	if(optimizeBaseFreq) {
		sum = .0;
		For(i,4) sum += FABS(unscaledBaseFreq[i]);
		For(i,4) baseFreq[i] = FABS(unscaledBaseFreq[i])/sum;
		
//		do {
//			sum = .0;
//			For(i,4) {
//				if(baseFreq[i] < 0.01) baseFreq[i]=0.01;
//				if(baseFreq[i] > 0.99) baseFreq[i]=0.99;
//				sum += baseFreq[i];
//			}
//			For(i,4) baseFreq[i]/=sum;
//		}
//		while((sum > 1.01) || (sum < 0.99));

//		do {
//			sum = .0;
//			For(i,4) {
//				if(baseFreq[i] < 0.001) baseFreq[i]=0.001;
//				if(baseFreq[i] > 0.999) baseFreq[i]=0.999;
//				sum += baseFreq[i];
//			}
//			For(i,4) baseFreq[i]/=sum;
//		}
//		while((sum > 1.001) || (sum < 0.999));

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

void Adjust_Min_Diff_Lk() {
	int exponent;
	exponent = (int)FLOOR(log10(FABS(lnlike)));
	if(sizeof(double) == 4) {
		min_diff_lk_local = POW(10.,exponent - FLT_DIG + 1);
    }
	if(sizeof(double) == 8) {
		min_diff_lk_local = POW(10.,exponent - DBL_DIG + 1);
    }
}

void CalcLnLike(double &lnlike, const std::vector<int> &weight, const std::vector<double> &condLike, const std::vector<double> &baseFreq, const std::vector<double> &pMatrix) {
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

#if POSTERIORMODE
#if DREAM
	lnlike += (log(6.0) - brlens/BRLENS_PR_MEAN - log(BRLENS_PR_MEAN));
#else
	lnlike += (log(6.0) - kappa - brlens/BRLENS_PR_MEAN - log(BRLENS_PR_MEAN));
#endif
#endif
}

double Lk() {
	Check_Br_Len_Bounds();
	NormalizeBaseFreq();
	CalcProbMat(pMatrix, baseFreq, kappa, brlens);
	CalcLnLike(lnlike, weight, condLike, baseFreq, pMatrix);
	Adjust_Min_Diff_Lk();
	return lnlike;
}

double Return_Abs_Lk() {
	Lk();
	return FABS(lnlike);
}

#define TOLX 1.0e-7
int Lnsrch(int n, double *xold, double fold, double *g, double *p, double *x, double *f, double stpmax, int *check, double (*func)()) {
	int i;
	double a, alam, alam2, alamin, b, disc, f2, fold2, rhs1, rhs2, slope, sum, temp, test, tmplam;
	double *local_xold;
	
	alam = alam2 = f2 = fold2 = tmplam = .0;
	
	local_xold = (double *)mCalloc(n,sizeof(double));
	For(i,n) local_xold[i] = xold[i];
	
	*check = 0;
	for(sum=0.0, i=0; i<n; i++)
		sum += p[i]*p[i];
	sum = SQRT(sum);
	if(sum > stpmax) {
		for(i=0; i<n; i++)
			p[i] *= stpmax/sum;
	}
	for(slope=0.0, i=0; i<n; i++)
		slope += g[i]*p[i];
	test = 0.0;
	for(i=0; i<n; i++) {
		temp = FABS(p[i])/MAX(FABS(local_xold[i]),1.0);
		if(temp > test)
			test = temp;
	}
	alamin = TOLX/test;
	alam = 1.0;
	for(;;) {
		for(i=0; i<n; i++) {
			x[i] = FABS(local_xold[i] + alam*p[i]);
			xold[i] = x[i];
		}
		if(i==n) {
			*f = (*func)();
		}
		else
			*f = 1. + fold + ALF*alam*slope;
		if(alam < alamin) {
			*check = 1;
			For(i,n) xold[i] = local_xold[i];
			Free(local_xold);
			return 0;
		}
		else if(*f <= fold+ALF*alam*slope) {
			For(i,n) xold[i] = local_xold[i];
			Free(local_xold);
			return 0;
		}
		else {
			if((alam < 1.0+SMALL) && (alam > 1.0-SMALL))
				tmplam = -slope/(2.0*(*f-fold-slope));
			else {
				rhs1 = *f-fold-alam*slope;
				rhs2 = f2-fold2-alam2*slope;
				a    = (rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				b    = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if(a < SMALL && a > -SMALL)
					tmplam = -slope/(2.0*b);
				else {
					disc = b*b-3.0*a*slope;
					if(disc<0.0)
						tmplam = 0.5*alam;
					else if(b <= 0.0)
						tmplam = (-b+SQRT(disc))/(3.0*a);
					else
						tmplam = -slope/(b+SQRT(disc));
				}
				if(tmplam>0.5*alam)
					tmplam = 0.5*alam;
			}
		}
		alam2 = alam;
		f2    = *f;
		fold2 = fold;
		alam  = MAX(tmplam,0.1*alam);
	}
	Free(local_xold);
	return 1;
}
#undef TOLX
#undef NRANSI

double Num_Derivatives_One_Param(double (*func)(), double f0, double *param, double stepsize, double *err, int precise) {
	int i, j;
	int n_iter;
	double errt, fac, hh, **a, ans;
	a = (double **)mCalloc(11,sizeof(double *));
	For(i,11) a[i] = (double *)mCalloc(11,sizeof(double));
	
	n_iter = 10;
	ans    = .0;
	
	if(stepsize < SMALL) {
		std::cerr << "\nh must be nonzero in Dfridr.\n";
		exit(1);
	}
	
	hh = stepsize;
	
	if(!precise) {
		*param   = *param+hh;
		a[0][0]  = (*func)();
		a[0][0] -= f0;
		a[0][0] /= hh;
		*param   = *param-hh;
		ans      = a[0][0];
	}
	else {
		*param   = *param+hh;
		a[0][0]  = (*func)();
		a[0][0] -= f0;
		a[0][0] /= hh;
		*param   = *param-hh;
		*err     = 1e30;
		for(i=1; i<n_iter; i++) {
			hh      /= 1.4;
			*param   = *param+hh;
			a[0][i]  = (*func)();
			a[0][i] -= f0;
			a[0][i] /= hh;
			*param   = *param-hh;
			fac      = 1.4*1.4;
			for(j=1; j<=i; j++) {
				a[j][i] = (a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
				fac     = 1.4*1.4*fac;
				errt    = MAX(FABS(a[j][i]-a[j-1][i]),FABS(a[j][i]-a[j-1][i-1]));
				if(errt <= *err) {
					*err = errt;
					ans  = a[j][i];
				}
			}
			if(FABS(a[i][i]-a[i-1][i-1]) >= 2.0*(*err))
				break;
		}
	}
	For(i,11) Free(a[i]);
	Free(a);
	
	return ans;
}

int Num_Derivative_Several_Param(double *param, int n_param, double stepsize, double (*func)(), double *derivatives) {
	int i;
	double err, f0;
	
	f0 = (*func)();
	
	For(i,n_param) {
		derivatives[i] = Num_Derivatives_One_Param(func, f0, param+i, stepsize, &err, 0);
	}
	
	return 1;
}

#define TOLX (4*EPS)
void BFGS(double *p,
		  int n,
		  double gtol,
		  double step_size,
		  double(*func)(),
		  int(*dfunc)(double *param,int n_param,double stepsize,double(*func)(),double *derivatives),
		  int(*lnsrch)(int n, double *xold, double fold, double *g, double *p, double *x, double *f, double stpmax, int *check, double(*func)()),
		  int *failed) {
	
	int check, i, its, j;
	double den, fac, fad, fae, fp, stpmax, sum=0.0, sumdg, sumxi, temp, test, fret;
	double *dg, *g, *hdg, **hessin, *pnew, *xi;
	
	hessin = (double **)mCalloc(n,sizeof(double *));
	For(i,n) hessin[i] = (double *)mCalloc(n,sizeof(double));
	dg   = (double *)mCalloc(n,sizeof(double ));
	g    = (double *)mCalloc(n,sizeof(double ));
	pnew = (double *)mCalloc(n,sizeof(double ));
	hdg  = (double *)mCalloc(n,sizeof(double ));
	xi   = (double *)mCalloc(n,sizeof(double ));
	
	fp = (*func)();
	(*dfunc)(p, n, step_size, func, g);
	
	for(i=0; i<n; i++) {
		for(j=0; j<n; j++)
			hessin[i][j] = 0.0;
		hessin[i][i] = 1.0;
		xi[i] = -g[i];
		sum  += p[i]*p[i];
	}
	
	stpmax = STPMX*MAX(SQRT(sum),(double)n);
	
	for(its=1; its<=ITMAX; its++) {
		lnsrch(n, p, fp, g, xi, pnew, &fret, stpmax, &check, func);
		
		fp = fret;
		
		for(i=0; i<n; i++) {
			xi[i] = pnew[i]-p[i];
			p[i]  = pnew[i];
		}
		
		test = 0.0;
		for(i=0; i<n; i++) {
			temp = FABS(xi[i])/MAX(FABS(p[i]),1.0);
			if(temp > test)
				test = temp;
		}
		if(test < TOLX) {
			(*func)();
			For(i,n) Free(hessin[i]);
			free(hessin);
			free(xi);
			free(pnew);
			free(hdg);
			free(g);
			free(dg);
			if(its == 1) {
				*failed = 1;
			}
			return;
		}
		
		for(i=0; i<n; i++)
			dg[i] = g[i];
		
		(*dfunc)(p, n, step_size, func, g);
		
		test = 0.0;
		den  = MAX(fret,1.0);
		for(i=0; i<n; i++) {
			temp = FABS(g[i])*MAX(FABS(p[i]),1.0)/den;
			if(temp > test)
				test = temp;
		}
		if(test < gtol) {
			(*func)();
			For(i,n) Free(hessin[i]);
			free(hessin);
			free(xi);
			free(pnew);
			free(hdg);
			free(g);
			free(dg);
			return;
		}
		
		for(i=0; i<n; i++)
			dg[i] = g[i]-dg[i];
		
		for(i=0; i<n; i++) {
			hdg[i] = 0.0;
			for(j=0; j<n; j++)
				hdg[i] += hessin[i][j]*dg[j];
		}
		
		fac = fae = sumdg = sumxi = 0.0;
		for(i=0; i<n; i++) {
			fac   += dg[i]*xi[i];
			fae   += dg[i]*hdg[i];
			sumdg += SQR(dg[i]);
			sumxi += SQR(xi[i]);
		}
		
		if(fac*fac > EPS*sumdg*sumxi) {
			fac = 1.0/fac;
			fad = 1.0/fae;
			for(i=0; i<n; i++)
				dg[i] = fac*xi[i]-fad*hdg[i];
			for(i=0; i<n; i++) {
				for(j=0; j<n; j++) {
					hessin[i][j] += fac*xi[i]*xi[j] - fad*hdg[i]*hdg[j] + fae*dg[i]*dg[j];
				}
			}
		}
		for(i=0; i<n; i++) {
			xi[i] = 0.0;
			for(j=0; j<n; j++)
				xi[i] -= hessin[i][j]*g[j];
		}
	}
	*failed = 1;
	For(i,n) Free(hessin[i]);
	free(hessin);
	free(xi);
	free(pnew);
	free(hdg);
	free(g);
	free(dg);
}
#undef TOLX

double Generic_Brent_Lk(double *param, double ax, double cx, double tol, int n_iter_max, int quickdirty, double(*obj_func)()) {
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
	fw = fv = fx = fu = -(*obj_func)();
	init_lnL = -fw;
	
	for(iter=1; iter<=BRENT_IT_MAX; iter++) {
		xm   = 0.5*(a+b);
		tol2 = 2.0*(tol1=tol*x+BRENT_ZEPS);
		
		if((fu > init_lnL + tol) && (quickdirty)) {
			(*param) = x;
			fu = (*obj_func)();
			return fu;
		}
		
		if((FABS(fu-old_lnL) < tol) || (iter > n_iter_max - 1)) {
			(*param) = x;
			fu = (*obj_func)();
			return fu;
		}
		
		if(FABS(e) > tol1) {
			r = (x-w)*(fx-fv);
			q = (x-v)*(fx-fw);
			p = (x-v)*q-(x-w)*r;
			q = 2.0*(q-r);
			if(q > 0.0)
				p = -p;
			q     = FABS(q);
			etemp = e;
			e     = d;
			if(FABS(p) >= FABS(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)) {
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
		
		u = (FABS(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		(*param) = u;
		old_lnL  = fu;
		fu       = -(*obj_func)();
		
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
			if(fu < fw || FABS(w-x) < SMALL) {
				v  = w;
				w  = u;
				fv = fw;
				fw = fu;
			}
			else if(fu < fv || FABS(v-x) < SMALL || FABS(v-w) < SMALL) {
				v  = u;
				fv = fu;
			}
		}
    }
	
	std::cerr << "\nToo many iterations in BRENT!\n";
	exit(1);
	return(-1);
}

void Optimize_State_Freqs() {
	optimizeBaseFreq = true;

	int i;
	For(i,4) {
		Generic_Brent_Lk(&(unscaledBaseFreq[i]),
						 0.,
						 100.,
						 min_diff_lk_local,
						 BRENT_IT_MAX,
						 0,
						 Lk);
	}
	
	optimizeBaseFreq = false;
}

void Optimize_TsTv() {
	Generic_Brent_Lk(&kappa,
					 0.001,
					 10.0,
					 min_diff_lk_local,
					 BRENT_IT_MAX,
					 0,
					 Lk);
}

void Optimiz_All_Free_Param() {
#if !DREAM
	Optimize_TsTv();
#endif
	Optimize_State_Freqs();
	Lk();
}

double Br_Len_Brent(double *brlens) {
	Generic_Brent_Lk(brlens,
					 brlensMin,
					 brlensMax,
					 min_diff_lk_local,
					 BRENT_IT_MAX,
					 0,
					 Lk);
	return lnlike;
}

void Optimize_Br_Len_Serie() {
	double lk_init = lnlike;
	Br_Len_Brent(&brlens);
	if(lnlike < lk_init - min_diff_lk_local) {
		std::cerr << "\nError in Optimize_Br_Len_Serie.\n";
		exit(1);
    }
}

void Round_Optimize(int n_round_max) {
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
		if(FABS(lk_new - lk_old) < min_diff_lk_local)
			break;
		else
			lk_old = lk_new;
		n_round++;
		each--;
    }
	
	Optimiz_All_Free_Param();
}

#endif
//================================================================================================= MuSE call
void CalcEvolDistance(std::string &s, std::vector<int> &q) {
	ResetParams();
	ReformatData(s, q, weight, condLike);
#if PHYML
	Round_Optimize(ROUND_MAX);
#endif
}

static inline void pileup_seq(const bam_pileup1_t *p, int pos, int ref_len, char **ref_ptr, std::string &s) {
	if(!p->is_del) {
		int c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos)];
		if(*ref_ptr) {
			int rb = pos < ref_len ? (*ref_ptr)[pos] : 'N';
			if(c=='=' || bam_nt16_table[c]==bam_nt16_table[rb]) {
				if(rb == 'N' || rb == 'n')
					s += '?';
				else
					s += char(toupper(rb));
			}
			else {
				if(c == 'N' || c == 'n')
					s += '?';
				else
					s += char(toupper(c));
			}
		}
		else {
			std::cout << "No reference.\n";
			exit(-1);
		}
	}
	else {
		s += '?';
	}
}

static int mplp_func(void *data, bam1_t *b) {
	extern uint8_t *bam_aux_get(const bam1_t *b, const char tag[2]);
	mplp_aux_t *ma = (mplp_aux_t*)data;
	int ret, skip = 0;
	do {
		ret = ma->iter ? bam_iter_read(ma->fp, ma->iter, b) : bam_read1(ma->fp, b);
		if(ret < 0)
			break;
		
		if(b->core.tid<0 || (b->core.flag&BAM_FUNMAP)) { // exclude unmapped reads
			skip = 1;
			continue;
		}

		int has_matched_ref = (*(ma->ref_id) == b->core.tid) ? 1 : 0;
		if(!has_matched_ref) {
			*(ma->ref_ptr)  = faidx_fetch_seq(ma->conf->fai, ma->h->target_name[b->core.tid], 0, 0x7fffffff, ma->ref_ptr, ma->ref_len);
			*(ma->ref_id)   = b->core.tid;
			has_matched_ref = 1;
		}
		
		if(ma->conf->flag & MPLP_ILLUMINA13) {
			uint8_t *qual = bam1_qual(b);
			for(int index = 0; index < b->core.l_qseq; ++index)
				qual[index] = qual[index] > 31 ? qual[index] - 31 : 0;
		}
		skip = 0;
	} while(skip);
	return ret;
}

static inline bool mySortFunc(std::pair<char,int> first, std::pair<char,int> second) {
	if(first.second > second.second)
		return true;
	else
		return false;
}

static int mpileup(mplp_conf_t *conf, int n, mplp_aux_t **data, bam_index_t *idxTumor, bam_index_t *idxNormal, char **ref_ptr, int &ref_tid, int &ref_len, std::ofstream &outf) {
	extern uint8_t *bam_aux_get(const bam1_t *b, const char tag[2]);
	extern int32_t bam_aux2i(const uint8_t *s);
	extern double kt_fisher_exact(int n11, int n12, int n21, int n22, double *_left, double *_right, double *two);

	int tid, pos, tid0 = -1, beg0 = 0, end0 = 1u<<29;

	int *n_plp;
	n_plp = (int*)calloc(n, sizeof(int*));

	const bam_pileup1_t **plp;
	plp = (const bam_pileup1_t **)calloc(n, sizeof(void*));

	bam_mplp_t iter;

	if(conf->reg) {
		int beg, end;		
		if(bam_parse_region(data[0]->h, conf->reg, &tid, &beg, &end) < 0) {
			fprintf(stderr, "[%s] malformatted region or wrong seqname for %d-th input.\n", __func__, 1);
			exit(1);
		}
		tid0 = tid, beg0 = beg, end0 = end;
		data[0]->iter = bam_iter_query(idxTumor, tid, beg, end);
		data[1]->iter = bam_iter_query(idxNormal, tid, beg, end);
	}

	if(tid0>=0 && tid0!=ref_tid) {
		*ref_ptr = faidx_fetch_seq(conf->fai, data[0]->h->target_name[tid0], 0, 0x7fffffff, ref_ptr, &ref_len);
		ref_tid = tid0;
	}

	iter = bam_mplp_init(n, mplp_func, (void**)data);
	bam_mplp_set_maxcnt(iter, conf->max_depth);

	while(bam_mplp_auto(iter, &tid, &pos, n_plp, plp) > 0) {
		// out of the region requested
		//
		if(conf->reg && (pos < beg0 || pos >= end0)) {
			continue;
		}
		
		// reload reference
		//
		if(tid != ref_tid) {
			*ref_ptr = faidx_fetch_seq(conf->fai, data[0]->h->target_name[tid], 0, 0x7fffffff, ref_ptr, &ref_len);
			ref_tid = tid;
		}
		
		refLetter = ((*ref_ptr) && pos < ref_len) ? toupper((*ref_ptr)[pos]) : '?';
		altLetter = ' ';
		
		// no need of calculation at the reference gap region
		//
		if(refLetter=='N' || refLetter=='n' || refLetter=='?') {
#if DEBUG
			std::cout << data[0]->h->target_name[tid] << '\t' << pos + 1 << '\t' << refLetter << '\t' << "reference gap region\n";
#endif
			continue;
		}

		// do not meet the minimum depth
		//
		if(n_plp[0] < conf->depth_cutoff || n_plp[1] < conf->depth_cutoff) {
#if DEBUG
			std::cout << data[0]->h->target_name[tid] << '\t' << pos + 1 << '\t' << refLetter << '\t' << "do not meet the minimum depth: n_plp[0]=" << n_plp[0] << " n_plp[1]=" << n_plp[1] << '\n';
#endif
			continue;
		}
		
		// raw pileup sequence
		//
		std::vector<std::string> rawS(n, "");
		std::vector< std::vector<int> > rawQ(n);
		rawQ[0].reserve(n_plp[0]);
		rawQ[1].reserve(n_plp[1]);
		
		// raw pileup element pointer
		//
		std::vector<const bam_pileup1_t*> rawTumorPileupElementPointer;
		rawTumorPileupElementPointer.reserve(n_plp[0]);
		std::vector<const bam_pileup1_t*> rawNormalPileupElementPointer;
		rawNormalPileupElementPointer.reserve(n_plp[1]);
		
		typedef std::map<std::string, pairPosition>::value_type valType;
		std::vector< std::map<std::string, pairPosition> > overlapPair(n);
		
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < n_plp[i]; j++) {
				const bam_pileup1_t *p = plp[i] + j;
				// generate pileup sequence
				//
				pileup_seq(p, pos, ref_len, ref_ptr, rawS[i]);
				int c_raw = bam1_qual(p->b)[p->qpos];
				if(c_raw > 93)
					rawQ[i].push_back(93);
				else
					rawQ[i].push_back(c_raw);

				if(i==0)
					rawTumorPileupElementPointer.push_back(p);
				else
					rawNormalPileupElementPointer.push_back(p);
				// overlap detection
				//
				std::string readName(bam1_qname(p->b));
				std::map<std::string, pairPosition>::iterator it = overlapPair[i].find(readName);
				if(it != overlapPair[i].end()) {
					(it->second).second = j;
				}
				else {
					pairPosition tmpPairPos;
					tmpPairPos.first = j, tmpPairPos.second = -1;
					overlapPair[i].insert(valType(readName, tmpPairPos));
				}
			}
		}

		// no alt allele in the raw pileup at all
		//
		if((std::count(rawS[0].begin(),rawS[0].end(),refLetter) + std::count(rawS[0].begin(),rawS[0].end(),'?')) == rawS[0].length()) {
				continue;
		}

		// remove overlap, follow mutect paper's suggestions
		//
		std::vector<const bam_pileup1_t*> noOverlapTumorPileupElementPointer(rawTumorPileupElementPointer);
		std::vector<const bam_pileup1_t*> noOverlapNormalPileupElementPointer(rawNormalPileupElementPointer);
		std::vector<std::string> tmpS(rawS);
		std::vector< std::vector<int> > tmpQ(rawQ);
		for(int i = 0; i < n; i++) {
			std::vector<int> removePos;
			for(std::map<std::string, pairPosition>::iterator it = overlapPair[i].begin(); it != overlapPair[i].end(); it++) {
				if((it->second).second != -1) {
					if(rawS[i][(it->second).first] != rawS[i][(it->second).second]) {
						if(i==0) {
							removePos.push_back((it->second).first);
							removePos.push_back((it->second).second);
						}
						else {
							if(rawS[i][(it->second).first] == refLetter)
								removePos.push_back((it->second).first);
							else if(rawS[i][(it->second).second] == refLetter)
								removePos.push_back((it->second).second);
						}
					}
					else {
						if(rawQ[i][(it->second).first] < rawQ[i][(it->second).second])
							removePos.push_back((it->second).first);
						else
							removePos.push_back((it->second).second);
					}
				}
			}
			if(!removePos.empty()) {
				std::sort(removePos.begin(), removePos.end());
				for(std::vector<int>::reverse_iterator it = removePos.rbegin(); it != removePos.rend(); it++) {
					tmpS[i].erase(tmpS[i].begin() + (*it));
					tmpQ[i].erase(tmpQ[i].begin() + (*it));
					if(i==0)
						noOverlapTumorPileupElementPointer.erase(noOverlapTumorPileupElementPointer.begin() + (*it));
					else
						noOverlapNormalPileupElementPointer.erase(noOverlapNormalPileupElementPointer.begin() + (*it));
				}
			}
		}

		// proximal gap filter
		//
		int  deletionCount  = 0;
		int  insertionCount = 0;
		for(int i = 0; i < noOverlapTumorPileupElementPointer.size(); i++) {
			// do the check so that there is no need of visiting all reads
			//
			if(deletionCount>=GAP_EVENT_CUTOFF || insertionCount>=GAP_EVENT_CUTOFF)
				break;
			// this is a deletion
			//
			if(noOverlapTumorPileupElementPointer[i]->is_del) {
				deletionCount++;
			}
			else {
				uint32_t *cigar = bam1_cigar(noOverlapTumorPileupElementPointer[i]->b);
				int positionInRead = noOverlapTumorPileupElementPointer[i]->qpos;
				int operationStart = 0;
				for(int j = 0; j < noOverlapTumorPileupElementPointer[i]->b->core.n_cigar; j++) {
					int operation = _cop(cigar[j]);
					int length    = _cln(cigar[j]);
					if(operation==BAM_CINS) {
						int distance = (operationStart>positionInRead) ? operationStart-positionInRead : positionInRead-(operationStart+length-1);
						if(distance<=GAP_EVENT_PROXIMITY) {
							insertionCount++;
							break;
						}
					}
					if(operation==BAM_CDEL) {
						int distance = (operationStart>positionInRead) ? operationStart-positionInRead : positionInRead-operationStart+1;
						if(distance<=GAP_EVENT_PROXIMITY) {
							deletionCount++;
							break;
						}
					}
					// move on to the next cigar operation
					//
					if(operation!=BAM_CDEL)
						operationStart += length;
				}
			}
		}
		if(deletionCount>=GAP_EVENT_CUTOFF || insertionCount>=GAP_EVENT_CUTOFF) {
#if DEBUG
			std::cout << data[0]->h->target_name[tid] << '\t' << pos + 1 << '\t' << refLetter << '\t' << "proximal gap\n";
#endif
			continue;
		}
		
		// key tumor pileup pointer, sequence and quality
		// and final normal pileup pointer
		// remove all '?', remove deletion and pass the MIN_QUALITY_SCORE, tumor data also remove MAPQ0 reads and mate unmapped reads
		//
		std::vector<std::string> s(n, "");
		std::vector< std::vector<int> > q(n);
		q[0].reserve(n_plp[0]);
		q[1].reserve(n_plp[1]);
		std::vector<const bam_pileup1_t*> noOverlapNoDelPassMinQNoMapQ0TumorPileupElementPointer;
		noOverlapNoDelPassMinQNoMapQ0TumorPileupElementPointer.reserve(n_plp[0]);
		for(int i = 0; i < noOverlapTumorPileupElementPointer.size(); i++) {
			const bam_pileup1_t *p = noOverlapTumorPileupElementPointer[i];
			if((tmpS[0][i]!='?') && (!(p->is_del)) && (p->b->core.flag&1) && (!(p->b->core.flag&8)) && (p->b->core.qual > 0) && (tmpQ[0][i] >= MIN_QUALITY_SCORE)) {
				noOverlapNoDelPassMinQNoMapQ0TumorPileupElementPointer.push_back(p);
				s[0] += tmpS[0][i];
				q[0].push_back(tmpQ[0][i]);
			}
		}
		std::vector<const bam_pileup1_t*> finalNormalPileupElementPointer;
		finalNormalPileupElementPointer.reserve(n_plp[1]);
		for(int i = 0; i < noOverlapNormalPileupElementPointer.size(); i++) {
			const bam_pileup1_t *p = noOverlapNormalPileupElementPointer[i];
			if((tmpS[1][i]!='?') && (!(p->is_del)) && (tmpQ[1][i] >= MIN_QUALITY_SCORE)) {
				finalNormalPileupElementPointer.push_back(p);
				s[1] += tmpS[1][i];
				q[1].push_back(tmpQ[1][i]);
			}
		}
	
		// recheck if the pileup meet the minimum depth
		//
		if(s[0].length() < conf->depth_cutoff || s[1].length() < conf->depth_cutoff) {
#if DEBUG
			std::cout << data[0]->h->target_name[tid] << '\t' << pos + 1 << '\t' << refLetter << '\t' << "do not meet the minimum depth\n";
#endif
			continue;
		}

		// identify what the somatic variant allele is
		//
		// used for alt_allele_in_normal post-filter
		//
		double obsTVAF = -1.0;
		double obsNVAF = -1.0;
		int nRefInTumor = 0;
		std::map<char, int> tumorVariantAlleleCount;
		tumorVariantAlleleCount.insert(std::map<char, int>::value_type('A', 0));
		tumorVariantAlleleCount.insert(std::map<char, int>::value_type('C', 0));
		tumorVariantAlleleCount.insert(std::map<char, int>::value_type('G', 0));
		tumorVariantAlleleCount.insert(std::map<char, int>::value_type('T', 0));

		for(int i = 0; i < s[0].length(); i++) {
			if(s[0][i]==refLetter) {
				nRefInTumor += 1;
			}
			else if(s[0][i]!='?') {
				tumorVariantAlleleCount[s[0][i]] += 1;
			}
		}
		
		// sort tumor variant allele count
		//
		std::vector< std::pair<char,int> > tumorVariantAlleleCountVec(tumorVariantAlleleCount.begin(), tumorVariantAlleleCount.end());
		std::sort(tumorVariantAlleleCountVec.begin(), tumorVariantAlleleCountVec.end(), &mySortFunc);

		// alt_allele_in_normal filter
		//
		if(tumorVariantAlleleCountVec[0].second != 0) {
			useSecondMutantAllele = false;
			// record the summary of normal data allele
			// key: ACGT
			// value: allele count
			// value: sum of quality scores
			//
			int nRefInNormal = 0;
			std::map<char, int> normalVariantAlleleCount;
			normalVariantAlleleCount.insert(std::map<char, int>::value_type('A', 0));
			normalVariantAlleleCount.insert(std::map<char, int>::value_type('C', 0));
			normalVariantAlleleCount.insert(std::map<char, int>::value_type('G', 0));
			normalVariantAlleleCount.insert(std::map<char, int>::value_type('T', 0));
			std::map<char, int> normalVariantAlleleQuality;
			normalVariantAlleleQuality.insert(std::map<char, int>::value_type('A', 0));
			normalVariantAlleleQuality.insert(std::map<char, int>::value_type('C', 0));
			normalVariantAlleleQuality.insert(std::map<char, int>::value_type('G', 0));
			normalVariantAlleleQuality.insert(std::map<char, int>::value_type('T', 0));
			for(int i = 0; i < s[1].length(); i++) {
				if(s[1][i]==refLetter) {
					nRefInNormal += 1;
				}
				else if(s[1][i]!='?') {
					normalVariantAlleleCount[s[1][i]]   += 1;
					normalVariantAlleleQuality[s[1][i]] += q[1][i];
				}
			}
			
			bool observedInNormal = false;

			// compare t_vaf with n_vaf
			//
			altLetter = tumorVariantAlleleCountVec[0].first;
			obsTVAF = double(tumorVariantAlleleCount[altLetter])/double(s[0].length());
			obsNVAF = double(normalVariantAlleleCount[altLetter])/double(s[1].length());

			// check if meet the rejection criteria
			//
			if((normalVariantAlleleCount[altLetter]>=NORMAL_OBSERVATION_COUNT || obsNVAF>=NORMAL_OBSERVATION_VAF) && normalVariantAlleleQuality[altLetter]>NORMAL_OBSERVATION_SUM_QUAL) {
				observedInNormal = true;
			}

			// check second largest variant allele if equal
			//
			if(observedInNormal && tumorVariantAlleleCountVec[1].second!=0) {
				char secondAllele = tumorVariantAlleleCountVec[1].first;
				obsTVAF = double(tumorVariantAlleleCount[secondAllele])/double(s[0].length());
				obsNVAF = double(normalVariantAlleleCount[secondAllele])/double(s[1].length());
				if(obsTVAF>=0.05) {
					if(!((normalVariantAlleleCount[secondAllele]>=NORMAL_OBSERVATION_COUNT || obsNVAF>=NORMAL_OBSERVATION_VAF) && normalVariantAlleleQuality[secondAllele]>NORMAL_OBSERVATION_SUM_QUAL)) {
						// need to make sure the reason for using the second somatic variant allele is because this is a germline variant position, either hetero- or homo-variant
						//
						// n11       n12         | null (hetero- or homo-variant)
						// n21       n22         | observation
						//-----------------------+-------------------------------
						// variant   non-variant | Total
						//
						double left, right, hetTwotail, homTwotail, prob;
						int coverageWithoutSecondAltAllele = s[1].length() - normalVariantAlleleCount[secondAllele];

						int contingencyTable[2][2] = {{0,0}, {0,0}};
						
						// test if this is a heterozygous variant position
						//
						contingencyTable[0][0] = coverageWithoutSecondAltAllele/2;
						contingencyTable[0][1] = coverageWithoutSecondAltAllele/2;
						contingencyTable[1][0] = normalVariantAlleleCount[altLetter];
						contingencyTable[1][1] = coverageWithoutSecondAltAllele - normalVariantAlleleCount[altLetter];
						prob = kt_fisher_exact(contingencyTable[0][0], contingencyTable[0][1], contingencyTable[1][0], contingencyTable[1][1], &left, &right, &hetTwotail);
						
						// test if this is a homozygous variant position
						//
						contingencyTable[0][0] = coverageWithoutSecondAltAllele;
						contingencyTable[0][1] = 0;
						contingencyTable[1][0] = normalVariantAlleleCount[altLetter];
						contingencyTable[1][1] = coverageWithoutSecondAltAllele - normalVariantAlleleCount[altLetter];
						prob = kt_fisher_exact(contingencyTable[0][0], contingencyTable[0][1], contingencyTable[1][0], contingencyTable[1][1], &left, &right, &homTwotail);

						if(hetTwotail>=GERMLINE_VARIANT_FET_CUTOFF || homTwotail>=GERMLINE_VARIANT_FET_CUTOFF) {
							observedInNormal = false;
							useSecondMutantAllele = true;
							altLetter = secondAllele;
						}
					}
				}
			}
			
			if(observedInNormal) {
#if DEBUG
				std::cout << data[0]->h->target_name[tid] << '\t' << pos + 1 << '\t' << refLetter << '\t' << "variant allele observed in normal data\n";
#endif
				continue;
			}
		}
		else {
			// no variant allele left in the tumor pileup after removing overlaps and generating the final pileup
			//
			continue;
		}

		// alternate allele frequency in the tumor less than the cutoff filter
		//
		if(obsTVAF < conf->minAltFraction) {
#if DEBUG
			std::cout << data[0]->h->target_name[tid] << '\t' << pos + 1 << '\t' << refLetter << '\t' << "observed variant allele fraction below cutoff\n";
#endif
			continue;
		}
		
		// strand bias post-filter by Fisher Exact Test
		//
		// use rawTumorPileupElementPointer and rawS to generate the contingency table
		// overlapped reads are kept but the criteria of forming the final pileup pointer are used, i.e., remove '?', deletion, less than MIN_QUALITY_SCORE, MAPQ0 and mate unmapped reads
		//
		// used for strand bias fisher exact test post-filter
		// n11     n12     | Alt
		// n21     n22     | Non-Alt
		//-----------------+-------
		// Forward Reverse | Total
		//
		int strandBiasContingencyTable[2][2] = {{0,0}, {0,0}};
		for(int i = 0; i < rawS[0].length(); i++) {
			const bam_pileup1_t *p = rawTumorPileupElementPointer[i];
			if((rawS[0][i]!='?') && (!(p->is_del)) && (p->b->core.flag&1) && (!(p->b->core.flag&8)) && (p->b->core.qual > 0) && (rawQ[0][i] >= MIN_QUALITY_SCORE)) {
				if(rawS[0][i]!=altLetter) {
					if(bam1_strand(p->b))
						strandBiasContingencyTable[1][1] += 1;
					else
						strandBiasContingencyTable[1][0] += 1;
				}
				else {
					if(bam1_strand(p->b))
						strandBiasContingencyTable[0][1] += 1;
					else
						strandBiasContingencyTable[0][0] += 1;
				}
			}
		}
		double left, right, twotail, prob;
		prob = kt_fisher_exact(strandBiasContingencyTable[0][0], strandBiasContingencyTable[0][1], strandBiasContingencyTable[1][0], strandBiasContingencyTable[1][1], &left, &right, &twotail);
		if(twotail <= STRAND_BIAS_FET_CUTOFF) {
#if DEBUG
			std::cout << data[0]->h->target_name[tid] << '\t' << pos + 1 << '\t' << refLetter << '\t' << "strand bias: twotail=" << twotail << '\n';
#endif
			continue;
		}

		// continue using the above key tumor pileup,
		// but now focus on the quality of reads supporting somatic variant allele, in order to see if it is possible to reject/filter this position
		// in the meanwhile, further filtration (similar to mutect but with modifications) is applied to the above key tumor pileup to get the final tumor pileup
		// 1) if has "XT" tag, cannot be "XT:A:M"
		// 2) < 30% of bases are soft-clipped
		// 3) sum of the quality scores of the mismatches <= 100
		// the difference is,
		// 4) does not include reads whose MAPQ<=10
		//
		// initialize and resetting some variables
		//
		std::vector<const bam_pileup1_t*> finalTumorPileupElementPointer;
		finalTumorPileupElementPointer.reserve(n_plp[0]);
		std::string tmpTumorS = s[0];
		s[0] = "";
		std::vector<int> tmpTumorQ(q[0]);
		q[0].clear();
		q[0].reserve(n_plp[0]);
		
		// record all tumor somatic variant reads' MAPQ
		// if the mean of all MAPQs <= 10, then this is not a good position
		//
		std::vector<int> tumorSomaticMutantReadMAPQ;
		tumorSomaticMutantReadMAPQ.reserve(n_plp[0]);
		
		// somatic variant clustered filter
		// this filter follows mutect setting that does not include reads that,
		// 1) if has "XT" tag, cannot be "XT:A:M"
		// 2) < 30% of bases are soft-clipped
		// 3) sum of the quality scores of the mismatches <= 100
		// the difference is,
		// 4) does not include reads whose MAPQ<=10
		// 5) the MEDIAN and MAD calculation is based on all appropriate reads' smallest distances
		//
		std::vector<int> smallestDistance;
		smallestDistance.reserve(n_plp[0]);

		// record the maximum somatic allele quality score and mapping score
		// the read must pass the following filtration:
		// 1) if has "XT" tag, cannot be "XT:A:M"
		// 2) < 30% of bases are soft-clipped
		// 3) sum of the quality scores of the mismatches <= 100
		// 4) if has "XT" tag, cannot be "XT:A:R"
		// 5) does not include reads whose MAPQ<=10
		// 6) read paired and read mapped in proper pair (mate mapped has been tested when generating the final pileup)
		//
		bool existRequiredRead = false;
		
		for(int i = 0; i < noOverlapNoDelPassMinQNoMapQ0TumorPileupElementPointer.size(); i++) {
			const bam_pileup1_t *p = noOverlapNoDelPassMinQNoMapQ0TumorPileupElementPointer[i];
			// mapping quality score
			//
			int mapQ = (int)(p->b->core.qual);
			if(tmpTumorS[i]==altLetter) {
				tumorSomaticMutantReadMAPQ.push_back(mapQ);
			}
			if(mapQ<=LOW_MAPQ)
				continue;
			// "XT:A:M" tag
			//
			uint8_t *symbol = bam_aux_get(p->b, "XT");
			if(symbol) {
				if(char(*(symbol+1))=='M')
					continue;
			}
			// < 30% of bases have been soft-clipped
			//
			uint32_t softclipLength = 0;
			uint32_t *cigar_for_sc  = bam1_cigar(p->b);
			for(int ithCigar = 0; ithCigar < p->b->core.n_cigar; ithCigar++) {
				int op = cigar_for_sc[ithCigar] & BAM_CIGAR_MASK;
				if(op == BAM_CSOFT_CLIP)
					softclipLength += cigar_for_sc[ithCigar] >> BAM_CIGAR_SHIFT;
			}
			if((double)softclipLength/(double)(p->b->core.l_qseq) >= HEAVY_SOFTCLIP)
				continue;
			// sum of the quality scores of the mismatches <= 100, here mismatches do not include insertion and deletion
			// use bam_fillmd1_core for reference
			//
			uint8_t     *seq   = bam1_seq(p->b);
			uint8_t     *qual  = bam1_qual(p->b);
			uint32_t    *cigar = bam1_cigar(p->b);
			bam1_core_t *core  = &(p->b->core);
			int k, x, y;
			int mismatchSumQual = 0;
			for(k=y=0, x = core->pos; k < core->n_cigar; ++k) {
				int j, l = _cln(cigar[k]), op = _cop(cigar[k]);
				if(op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
					for(j = 0; j < l; ++j) {
						int z = y + j;
						int c1 = bam1_seqi(seq, z), c2 = bam_nt16_table[(int)(*ref_ptr)[x+j]];
						if((*ref_ptr)[x+j] == 0)
							break; // out of boundary
						if((c1 == c2 && c1 != 15 && c2 != 15) || c1 == 0) { // a match
							continue;
						}
						else {
							mismatchSumQual += int(qual[z]);
						}
					}
					if(j < l)
						break;
					x += l;
					y += l;
				}
				else if(op == BAM_CDEL) {
					for(j = 0; j < l; ++j) {
						if((*ref_ptr)[x+j] == 0)
							break;
					}
					if(j < l)
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
			// now the read passes above filters
			//
			finalTumorPileupElementPointer.push_back(p);
			s[0] += tmpTumorS[i];
			q[0].push_back(tmpTumorQ[i]);
			// collect smallest distance if it is a somatic variant
			//
			if(tmpTumorS[i]==altLetter) {
				int rightmostCoordinate = bam_calend(&(p->b->core), bam1_cigar(p->b));
				int toRightmostDistance = abs(pos - rightmostCoordinate);
				int leftmostCoordinate  = p->b->core.pos;
				int toLeftmostDistance  = abs(pos - leftmostCoordinate);
				int minDistance = (toRightmostDistance >= toLeftmostDistance ? toLeftmostDistance : toRightmostDistance);
				smallestDistance.push_back(minDistance);
				// two more tests to decide if record the read MAPQ and base quality
				// "XT:A:R" tag
				//
				if(symbol) {
					if(char(*(symbol+1))=='R')
						continue;
				}
				// read paired and read mapped in proper pair
				//
				if((p->b->core.flag&1) && (p->b->core.flag&2)) {
					if(mapQ >= REQUIRED_MAPQ && tmpTumorQ[i] >= REQUIRED_BASE_QUALITY) {
						existRequiredRead = true;
					}
				}
			}
		}

		// not meet the required MAPQ or base quality
		//
		if(!existRequiredRead) {
#if DEBUG
			std::cout << data[0]->h->target_name[tid] << '\t' << pos + 1 << '\t' << refLetter << '\t' << "no read meets the required MAPQ and base quality\n";
#endif
			continue;
		}
		
		// recheck if there is no alt allele in the final tumor pileup
		//
		if(std::count(s[0].begin(),s[0].end(),altLetter) == 0) {
#if DEBUG
			std::cout << data[0]->h->target_name[tid] << '\t' << pos + 1 << '\t' << refLetter << '\t' << "no appropriate somatic variant allele\n";
#endif
			continue;
		}
		
		// low_mean_mapq filter
		// calculate mean
		//
		double sumMAPQ = 0.0;
		for(int i = 0; i < tumorSomaticMutantReadMAPQ.size(); i++) {
			sumMAPQ += (double)tumorSomaticMutantReadMAPQ[i];
		}
		double meanMAPQ = sumMAPQ/(double)(tumorSomaticMutantReadMAPQ.size());
		if(meanMAPQ <= LOW_MEAN_MAPQ) {
#if DEBUG
			std::cout << data[0]->h->target_name[tid] << '\t' << pos + 1 << '\t' << refLetter << '\t' << "low average read mapq\n";
#endif
			continue;
		}
		
		// clustered position post-filter
		// calculate median
		//
		int nMutant = (int)(smallestDistance.size());
		std::sort(smallestDistance.begin(), smallestDistance.end());
		double distanceMedian = -1.0;
		if(nMutant%2 == 1) {
			distanceMedian = (double)smallestDistance[(int)floor(nMutant/2)];
		}
		else {
			double distanceLowerMiddle = (double)smallestDistance[(int)(nMutant/2)];
			double distanceUpperMiddle = (double)smallestDistance[(int)(nMutant/2) - 1];
			distanceMedian = (distanceLowerMiddle + distanceUpperMiddle)/2.0;
		}
		// calculate median absolute deviation (MAD)
		//
		std::vector<double> distanceDev(nMutant, 0.0);
		for(int i = 0; i < nMutant; i++) {
			distanceDev[i] = fabs((double)smallestDistance[i] - distanceMedian);
		}
		std::sort(distanceDev.begin(), distanceDev.end());
		double distanceMAD = -1.0;
		if(nMutant%2 == 1) {
			distanceMAD = distanceDev[(int)floor(nMutant/2)];
		}
		else {
			double distanceDevLowerMiddle = distanceDev[(int)(nMutant/2)];
			double distanceDevUpperMiddle = distanceDev[(int)(nMutant/2) - 1];
			distanceMAD = (distanceDevLowerMiddle + distanceDevUpperMiddle)/2.0;
		}

		if(distanceMedian<=CLUSTERED_POSITION_MEDIAN_CUTOFF && distanceMAD<=CLUSTERED_POSITION_MAD_CUTOFF) {
#if DEBUG
			std::cout << data[0]->h->target_name[tid] << '\t' << pos + 1 << '\t' << refLetter << '\t' << "clustered position\n";
#endif
			continue;
		}

		// if a position passes all the above quality control tests, then calculate the evolutionary distance
		// before doing the caculation, collect read group information from the final tumor pileup
		//
		// has some problems with char*, so convert it into std:string
		//
		std::set<std::string> ReadGroupInfo;
		for(int i = 0; i < s[0].length(); i++) {
			if(s[0][i]==altLetter) {
				const uint8_t *rg = bam_aux_get(finalTumorPileupElementPointer[i]->b, "RG");
				if(rg) {
					std::string rgContent((char*)(rg+1));
					if(std::find(ReadGroupInfo.begin(), ReadGroupInfo.end(), rgContent)==ReadGroupInfo.end()) {
						ReadGroupInfo.insert(rgContent);
					}
				}
			}
		}
		
		// count the number of reliable reference allele in the final normal pileup in order to overcome the possibility of normal data missing germ-line variant because of the low normal data coverage
		// the cutoff will be chosen at 10 reference alleles, but could be more stringent, like 20
		//
		int nGoodRefInNormal = 0;
		for(int i = 0; i < s[1].length(); i++) {
			if(s[1][i]==refLetter) {
				// Mapping Quality score > 0
				//
				if(finalNormalPileupElementPointer[i]->b->core.qual > 0) {
					// read not mapped in proper pair
					//
					if((finalNormalPileupElementPointer[i]->b->core.flag&1) && (finalNormalPileupElementPointer[i]->b->core.flag&2)) {
						// mate not unmapped
						//
						if(!(finalNormalPileupElementPointer[i]->b->core.flag&8))
							nGoodRefInNormal += 1;
					}
				}
			}
		}
		
		CalcEvolDistance(s[0], q[0]);

		if(brlens > conf->min_output_brlens) {
			std::vector<double> tumorBaseFreq(baseFreq);
			double tumorKappa  = kappa;
			double tumorBrlens = brlens;
			if(std::count(s[1].begin(), s[1].end(), altLetter) != 0) {
				CalcEvolDistance(s[1], q[1]);
				outf << data[0]->h->target_name[tid] << '\t' << pos + 1 << '\t' << refLetter << '\t' << altLetter << '\t' << s[0].length() << '\t' << s[1].length() << '\t';
				outf << std::scientific << obsTVAF << '\t' << tumorBaseFreq[0] << '\t' << tumorBaseFreq[1] << '\t' << tumorBaseFreq[2] << '\t' << tumorBaseFreq[3] << '\t' << tumorKappa << '\t' << tumorBrlens << '\t';
				outf << std::scientific << obsNVAF << '\t' << baseFreq[0] << '\t' << baseFreq[1] << '\t' << baseFreq[2] << '\t' << baseFreq[3] << '\t' << kappa << '\t' << brlens << '\t' << nGoodRefInNormal << '\t' << ReadGroupInfo.size() << '\t';
				if(useSecondMutantAllele)
					outf << "Y" << '\t';
				else
					outf << "N" << '\t';
				outf << std::count(s[0].begin(), s[0].end(), refLetter) << '\t' << std::count(s[0].begin(), s[0].end(), altLetter) << '\t' << std::count(s[1].begin(), s[1].end(), refLetter) << '\t' << std::count(s[1].begin(), s[1].end(), altLetter) << '\t';
			}
			else {
				outf << data[0]->h->target_name[tid] << '\t' << pos + 1 << '\t' << refLetter << '\t' << altLetter << '\t' << s[0].length() << '\t' << s[1].length() << '\t';
				outf << std::scientific << obsTVAF << '\t' << tumorBaseFreq[0] << '\t' << tumorBaseFreq[1] << '\t' << tumorBaseFreq[2] << '\t' << tumorBaseFreq[3] << '\t' << tumorKappa << '\t' << tumorBrlens << '\t';
				outf << std::scientific << obsNVAF << '\t';
				outf << "NA" << '\t' << "NA" << '\t' << "NA" << '\t' << "NA" << '\t' << "NA" << '\t' << MIN_BRLENS << '\t' << nGoodRefInNormal << '\t' << ReadGroupInfo.size() << '\t';
				if(useSecondMutantAllele)
					outf << "Y" << '\t';
				else
					outf << "N" << '\t';
				outf << std::count(s[0].begin(), s[0].end(), refLetter) << '\t' << std::count(s[0].begin(), s[0].end(), altLetter) << '\t' << std::count(s[1].begin(), s[1].end(), refLetter) << '\t' << std::count(s[1].begin(), s[1].end(), altLetter) << '\t';
			}
			// below is added content according to the requirements of TCGA VCF Specification
			//
			std::string altString = "";
			std::string normalFormat = "";
			std::string tumorFormat  = "";

			altString += altLetter;
			
			// normal raw data
			//
			std::map<char, int> rawNormalAlleleCount;
			rawNormalAlleleCount.insert(std::map<char, int>::value_type('A', 0));
			rawNormalAlleleCount.insert(std::map<char, int>::value_type('C', 0));
			rawNormalAlleleCount.insert(std::map<char, int>::value_type('G', 0));
			rawNormalAlleleCount.insert(std::map<char, int>::value_type('T', 0));
			rawNormalAlleleCount.insert(std::map<char, int>::value_type('?', 0));
			std::map<char, int> rawNormalAlleleQuality;
			rawNormalAlleleQuality.insert(std::map<char, int>::value_type('A', 0));
			rawNormalAlleleQuality.insert(std::map<char, int>::value_type('C', 0));
			rawNormalAlleleQuality.insert(std::map<char, int>::value_type('G', 0));
			rawNormalAlleleQuality.insert(std::map<char, int>::value_type('T', 0));
			rawNormalAlleleQuality.insert(std::map<char, int>::value_type('?', 0));
			for(int i = 0; i < rawS[1].length(); i++) {
				rawNormalAlleleCount[rawS[1][i]]   += 1;
				rawNormalAlleleQuality[rawS[1][i]] += rawQ[1][i];
			}

			// tumor raw data
			//
			std::map<char, int> rawTumorAlleleCount;
			rawTumorAlleleCount.insert(std::map<char, int>::value_type('A', 0));
			rawTumorAlleleCount.insert(std::map<char, int>::value_type('C', 0));
			rawTumorAlleleCount.insert(std::map<char, int>::value_type('G', 0));
			rawTumorAlleleCount.insert(std::map<char, int>::value_type('T', 0));
			rawTumorAlleleCount.insert(std::map<char, int>::value_type('?', 0));
			std::map<char, int> rawTumorAlleleQuality;
			rawTumorAlleleQuality.insert(std::map<char, int>::value_type('A', 0));
			rawTumorAlleleQuality.insert(std::map<char, int>::value_type('C', 0));
			rawTumorAlleleQuality.insert(std::map<char, int>::value_type('G', 0));
			rawTumorAlleleQuality.insert(std::map<char, int>::value_type('T', 0));
			rawTumorAlleleQuality.insert(std::map<char, int>::value_type('?', 0));
			for(int i = 0; i < rawS[0].length(); i++) {
				rawTumorAlleleCount[rawS[0][i]]   += 1;
				rawTumorAlleleQuality[rawS[0][i]] += rawQ[0][i];
			}
			
			if((rawNormalAlleleCount[refLetter]+rawNormalAlleleCount[altLetter]+rawNormalAlleleCount['?']) == rawS[1].length()) {
				// it is not a germline position based on the normal data
				//
				// normal Format first
				//
				normalFormat = "0/0:" + IntToString(rawS[1].length()-rawNormalAlleleCount['?']) + ":" + IntToString(rawNormalAlleleCount[refLetter]) + "," + IntToString(rawNormalAlleleCount[altLetter]) + ":";
				if(rawNormalAlleleCount[refLetter]==0) {
					normalFormat += "0,";
				}
				else {
					int tmp = round((double)rawNormalAlleleQuality[refLetter]/(double)rawNormalAlleleCount[refLetter]);
					normalFormat += (IntToString(tmp) + ",");
				}
				if(rawNormalAlleleCount[altLetter]==0) {
					normalFormat += "0:.";
				}
				else {
					int tmp = round((double)rawNormalAlleleQuality[altLetter]/(double)rawNormalAlleleCount[altLetter]);
					normalFormat += (IntToString(tmp) + ":.");
				}
				// tumor Format
				// always assume the GT is 0/1 unless test result shows otherwise
				//
				// n11       n12         | null (homo-variant)
				// n21       n22         | observation
				//-----------------------+-------------------------------
				// variant   ref         | Total
				//
				double left, right, homTwotail, prob;
				int coverage = rawTumorAlleleCount[refLetter] + rawTumorAlleleCount[altLetter];
				int contingencyTable[2][2] = {{0,0}, {0,0}};
				// test if this is a homozygous variant position
				//
				contingencyTable[0][0] = coverage;
				contingencyTable[0][1] = 0;
				contingencyTable[1][0] = rawTumorAlleleCount[altLetter];
				contingencyTable[1][1] = rawTumorAlleleCount[refLetter];
				prob = kt_fisher_exact(contingencyTable[0][0], contingencyTable[0][1], contingencyTable[1][0], contingencyTable[1][1], &left, &right, &homTwotail);
				
				if(homTwotail>=GERMLINE_VARIANT_FET_CUTOFF) {
					// the test result is not significant, so I cannot reject the hypothesis that this position is homo-variant
					//
					tumorFormat = "1/1:";
				}
				else {
					tumorFormat = "0/1:";
				}
				tumorFormat += (IntToString(rawS[0].length()-rawTumorAlleleCount['?']) + ":" + IntToString(rawTumorAlleleCount[refLetter]) + "," + IntToString(rawTumorAlleleCount[altLetter]) + ":");

				if(rawTumorAlleleCount[refLetter]==0) {
					tumorFormat += "0,";
				}
				else {
					int tmp = round((double)rawTumorAlleleQuality[refLetter]/(double)rawTumorAlleleCount[refLetter]);
					tumorFormat += (IntToString(tmp) + ",");
				}
				if(rawTumorAlleleCount[altLetter]==0) {
					tumorFormat += "0:2";
				}
				else {
					int tmp = round((double)rawTumorAlleleQuality[altLetter]/(double)rawTumorAlleleCount[altLetter]);
					tumorFormat += (IntToString(tmp) + ":2");
				}
			}
			else {
				// it might be a germline position, but it must meet the following requirements
				// 1) there is an allele (not refLetter or altLetter) count that is larger than altLetter's count, only check the largest count from those of the two not refLetter or altLetter
				// 2) fisher exact test can decide this is a heter-variant or homo-variant position
				//
				std::vector< std::pair<char,int> > rawNormalAlleleCountVec(rawNormalAlleleCount.begin(), rawNormalAlleleCount.end());
				std::sort(rawNormalAlleleCountVec.begin(), rawNormalAlleleCountVec.end(), &mySortFunc);
				char potentialGermlineAllele = ' ';
				for(int i = 0; i < rawNormalAlleleCountVec.size(); i++) {
					if(rawNormalAlleleCountVec[i].first!=refLetter && rawNormalAlleleCountVec[i].first!=altLetter && rawNormalAlleleCountVec[i].first!='?' && rawNormalAlleleCountVec[i].second>rawNormalAlleleCount[altLetter]) {
						potentialGermlineAllele = rawNormalAlleleCountVec[i].first;
						break;
					}
				}
				if(potentialGermlineAllele==' ') {
					// it is not a germline position based on the normal data
					//
					// normal Format first
					//
					normalFormat = "0/0:" + IntToString(rawS[1].length()-rawNormalAlleleCount['?']) + ":" + IntToString(rawNormalAlleleCount[refLetter]) + "," + IntToString(rawNormalAlleleCount[altLetter]) + ":";
					if(rawNormalAlleleCount[refLetter]==0) {
						normalFormat += "0,";
					}
					else {
						int tmp = round((double)rawNormalAlleleQuality[refLetter]/(double)rawNormalAlleleCount[refLetter]);
						normalFormat += (IntToString(tmp) + ",");
					}
					if(rawNormalAlleleCount[altLetter]==0) {
						normalFormat += "0:.";
					}
					else {
						int tmp = round((double)rawNormalAlleleQuality[altLetter]/(double)rawNormalAlleleCount[altLetter]);
						normalFormat += (IntToString(tmp) + ":.");
					}
					// tumor Format
					// always assume the GT is 0/1 unless test result shows otherwise
					//
					// n11       n12         | null (homo-variant)
					// n21       n22         | observation
					//-----------------------+-------------------------------
					// variant   ref         | Total
					//
					double left, right, homTwotail, prob;
					int coverage = rawTumorAlleleCount[refLetter] + rawTumorAlleleCount[altLetter];
					int contingencyTable[2][2] = {{0,0}, {0,0}};
					// test if this is a homozygous variant position
					//
					contingencyTable[0][0] = coverage;
					contingencyTable[0][1] = 0;
					contingencyTable[1][0] = rawTumorAlleleCount[altLetter];
					contingencyTable[1][1] = rawTumorAlleleCount[refLetter];
					prob = kt_fisher_exact(contingencyTable[0][0], contingencyTable[0][1], contingencyTable[1][0], contingencyTable[1][1], &left, &right, &homTwotail);
					
					if(homTwotail>=GERMLINE_VARIANT_FET_CUTOFF) {
						// the test result is not significant, so I cannot reject the hypothesis that this position is homo-variant
						//
						tumorFormat = "1/1:";
					}
					else {
						tumorFormat = "0/1:";
					}
					tumorFormat += (IntToString(rawS[0].length()-rawTumorAlleleCount['?']) + ":" + IntToString(rawTumorAlleleCount[refLetter]) + "," + IntToString(rawTumorAlleleCount[altLetter]) + ":");
					
					if(rawTumorAlleleCount[refLetter]==0) {
						tumorFormat += "0,";
					}
					else {
						int tmp = round((double)rawTumorAlleleQuality[refLetter]/(double)rawTumorAlleleCount[refLetter]);
						tumorFormat += (IntToString(tmp) + ",");
					}
					if(rawTumorAlleleCount[altLetter]==0) {
						tumorFormat += "0:2";
					}
					else {
						int tmp = round((double)rawTumorAlleleQuality[altLetter]/(double)rawTumorAlleleCount[altLetter]);
						tumorFormat += (IntToString(tmp) + ":2");
					}
				}
				else {
					// n11       n12         | null (hetero- or homo-variant)
					// n21       n22         | observation
					//-----------------------+-------------------------------
					// variant   ref         | Total
					//
					double left, right, hetTwotail, homTwotail, prob;
					int coverage = rawNormalAlleleCount[refLetter] + rawNormalAlleleCount[potentialGermlineAllele];
					int contingencyTable[2][2] = {{0,0}, {0,0}};
					
					// test if this is a heterozygous variant position
					//
					contingencyTable[0][0] = coverage/2;
					contingencyTable[0][1] = coverage/2;
					contingencyTable[1][0] = rawNormalAlleleCount[potentialGermlineAllele];
					contingencyTable[1][1] = rawNormalAlleleCount[refLetter];
					prob = kt_fisher_exact(contingencyTable[0][0], contingencyTable[0][1], contingencyTable[1][0], contingencyTable[1][1], &left, &right, &hetTwotail);
					
					// test if this is a homozygous variant position
					//
					contingencyTable[0][0] = coverage;
					contingencyTable[0][1] = 0;
					contingencyTable[1][0] = rawNormalAlleleCount[potentialGermlineAllele];
					contingencyTable[1][1] = rawNormalAlleleCount[refLetter];
					prob = kt_fisher_exact(contingencyTable[0][0], contingencyTable[0][1], contingencyTable[1][0], contingencyTable[1][1], &left, &right, &homTwotail);
					
					if(hetTwotail>=GERMLINE_VARIANT_FET_CUTOFF) {
						std::string tmpS(1, potentialGermlineAllele);
						altString += ("," + tmpS);
						
						// normal Format first
						//
						normalFormat = "0/2:" + IntToString(rawS[1].length()-rawNormalAlleleCount['?']) + ":" + IntToString(rawNormalAlleleCount[refLetter]) + "," + IntToString(rawNormalAlleleCount[altLetter]) + "," + IntToString(rawNormalAlleleCount[potentialGermlineAllele]) + ":";
						if(rawNormalAlleleCount[refLetter]==0) {
							normalFormat += "0,";
						}
						else {
							int tmp = round((double)rawNormalAlleleQuality[refLetter]/(double)rawNormalAlleleCount[refLetter]);
							normalFormat += (IntToString(tmp) + ",");
						}
						if(rawNormalAlleleCount[altLetter]==0) {
							normalFormat += "0,";
						}
						else {
							int tmp = round((double)rawNormalAlleleQuality[altLetter]/(double)rawNormalAlleleCount[altLetter]);
							normalFormat += (IntToString(tmp) + ",");
						}
						if(rawNormalAlleleCount[potentialGermlineAllele]==0) {
							normalFormat += "0:.";
						}
						else {
							int tmp = round((double)rawNormalAlleleQuality[potentialGermlineAllele]/(double)rawNormalAlleleCount[potentialGermlineAllele]);
							normalFormat += (IntToString(tmp) + ":.");
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
						double left, right, homTwotail, prob;
						int coverage = rawTumorAlleleCount[refLetter] + rawTumorAlleleCount[altLetter];
						int contingencyTable[2][2] = {{0,0}, {0,0}};
						// test if this is a homozygous variant position
						//
						contingencyTable[0][0] = coverage;
						contingencyTable[0][1] = 0;
						contingencyTable[1][0] = rawTumorAlleleCount[altLetter];
						contingencyTable[1][1] = rawTumorAlleleCount[refLetter];
						prob = kt_fisher_exact(contingencyTable[0][0], contingencyTable[0][1], contingencyTable[1][0], contingencyTable[1][1], &left, &right, &homTwotail);
						
						if(homTwotail>=GERMLINE_VARIANT_FET_CUTOFF) {
							// the test result is not significant, so I cannot reject the hypothesis that this position is homo-variant
							//
							if(rawTumorAlleleCount[potentialGermlineAllele]!=0) {
								tumorFormat = "1/1/2:";
							}
							else {
								tumorFormat = "1/1:";
							}
						}
						else {
							if(rawTumorAlleleCount[potentialGermlineAllele]!=0) {
								tumorFormat = "0/1/2:";
							}
							else {
								tumorFormat = "0/1:";
							}
						}
						tumorFormat += (IntToString(rawS[0].length()-rawTumorAlleleCount['?']) + ":" + IntToString(rawTumorAlleleCount[refLetter]) + "," + IntToString(rawTumorAlleleCount[altLetter]) + "," + IntToString(rawTumorAlleleCount[potentialGermlineAllele]) + ":");
						
						if(rawTumorAlleleCount[refLetter]==0) {
							tumorFormat += "0,";
						}
						else {
							int tmp = round((double)rawTumorAlleleQuality[refLetter]/(double)rawTumorAlleleCount[refLetter]);
							tumorFormat += (IntToString(tmp) + ",");
						}
						if(rawTumorAlleleCount[altLetter]==0) {
							tumorFormat += "0,";
						}
						else {
							int tmp = round((double)rawTumorAlleleQuality[altLetter]/(double)rawTumorAlleleCount[altLetter]);
							tumorFormat += (IntToString(tmp) + ",");
						}
						if(rawTumorAlleleCount[potentialGermlineAllele]==0) {
							tumorFormat += "0:2";
						}
						else {
							int tmp = round((double)rawTumorAlleleQuality[potentialGermlineAllele]/(double)rawTumorAlleleCount[potentialGermlineAllele]);
							tumorFormat += (IntToString(tmp) + ":2");
						}
					}
					else if(homTwotail>=GERMLINE_VARIANT_FET_CUTOFF) {
						std::string tmpS(1, potentialGermlineAllele);
						altString += ("," + tmpS);
						
						// normal Format first
						//
						normalFormat = "2/2:" + IntToString(rawS[1].length()-rawNormalAlleleCount['?']) + ":" + IntToString(rawNormalAlleleCount[refLetter]) + "," + IntToString(rawNormalAlleleCount[altLetter]) + "," + IntToString(rawNormalAlleleCount[potentialGermlineAllele]) + ":";
						if(rawNormalAlleleCount[refLetter]==0) {
							normalFormat += "0,";
						}
						else {
							int tmp = round((double)rawNormalAlleleQuality[refLetter]/(double)rawNormalAlleleCount[refLetter]);
							normalFormat += (IntToString(tmp) + ",");
						}
						if(rawNormalAlleleCount[altLetter]==0) {
							normalFormat += "0,";
						}
						else {
							int tmp = round((double)rawNormalAlleleQuality[altLetter]/(double)rawNormalAlleleCount[altLetter]);
							normalFormat += (IntToString(tmp) + ",");
						}
						if(rawNormalAlleleCount[potentialGermlineAllele]==0) {
							normalFormat += "0:.";
						}
						else {
							int tmp = round((double)rawNormalAlleleQuality[potentialGermlineAllele]/(double)rawNormalAlleleCount[potentialGermlineAllele]);
							normalFormat += (IntToString(tmp) + ":.");
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
						double left, right, homTwotail, prob;
						int coverage = rawTumorAlleleCount[refLetter] + rawTumorAlleleCount[altLetter];
						int contingencyTable[2][2] = {{0,0}, {0,0}};
						// test if this is a homozygous variant position
						//
						contingencyTable[0][0] = coverage;
						contingencyTable[0][1] = 0;
						contingencyTable[1][0] = rawTumorAlleleCount[altLetter];
						contingencyTable[1][1] = rawTumorAlleleCount[refLetter];
						prob = kt_fisher_exact(contingencyTable[0][0], contingencyTable[0][1], contingencyTable[1][0], contingencyTable[1][1], &left, &right, &homTwotail);
						
						if(homTwotail>=GERMLINE_VARIANT_FET_CUTOFF) {
							// the test result is not significant, so I cannot reject the hypothesis that this position is homo-variant
							//
							if(rawTumorAlleleCount[potentialGermlineAllele]!=0) {
								tumorFormat = "1/1/2:";
							}
							else {
								tumorFormat = "1/1:";
							}
						}
						else {
							if(rawTumorAlleleCount[potentialGermlineAllele]!=0) {
								tumorFormat = "0/1/2:";
							}
							else {
								tumorFormat = "0/1:";
							}
						}
						tumorFormat += (IntToString(rawS[0].length()-rawTumorAlleleCount['?']) + ":" + IntToString(rawTumorAlleleCount[refLetter]) + "," + IntToString(rawTumorAlleleCount[altLetter]) + "," + IntToString(rawTumorAlleleCount[potentialGermlineAllele]) + ":");
						
						if(rawTumorAlleleCount[refLetter]==0) {
							tumorFormat += "0,";
						}
						else {
							int tmp = round((double)rawTumorAlleleQuality[refLetter]/(double)rawTumorAlleleCount[refLetter]);
							tumorFormat += (IntToString(tmp) + ",");
						}
						if(rawTumorAlleleCount[altLetter]==0) {
							tumorFormat += "0,";
						}
						else {
							int tmp = round((double)rawTumorAlleleQuality[altLetter]/(double)rawTumorAlleleCount[altLetter]);
							tumorFormat += (IntToString(tmp) + ",");
						}
						if(rawTumorAlleleCount[potentialGermlineAllele]==0) {
							tumorFormat += "0:2";
						}
						else {
							int tmp = round((double)rawTumorAlleleQuality[potentialGermlineAllele]/(double)rawTumorAlleleCount[potentialGermlineAllele]);
							tumorFormat += (IntToString(tmp) + ":2");
						}
					}
					else {
						// it is not a germline position based on the normal data
						//
						// normal Format first
						//
						normalFormat = "0/0:" + IntToString(rawS[1].length()-rawNormalAlleleCount['?']) + ":" + IntToString(rawNormalAlleleCount[refLetter]) + "," + IntToString(rawNormalAlleleCount[altLetter]) + ":";
						if(rawNormalAlleleCount[refLetter]==0) {
							normalFormat += "0,";
						}
						else {
							int tmp = round((double)rawNormalAlleleQuality[refLetter]/(double)rawNormalAlleleCount[refLetter]);
							normalFormat += (IntToString(tmp) + ",");
						}
						if(rawNormalAlleleCount[altLetter]==0) {
							normalFormat += "0:.";
						}
						else {
							int tmp = round((double)rawNormalAlleleQuality[altLetter]/(double)rawNormalAlleleCount[altLetter]);
							normalFormat += (IntToString(tmp) + ":.");
						}
						// tumor Format
						// always assume the GT is 0/1 unless test result shows otherwise
						//
						// n11       n12         | null (homo-variant)
						// n21       n22         | observation
						//-----------------------+-------------------------------
						// variant   ref         | Total
						//
						double left, right, homTwotail, prob;
						int coverage = rawTumorAlleleCount[refLetter] + rawTumorAlleleCount[altLetter];
						int contingencyTable[2][2] = {{0,0}, {0,0}};
						// test if this is a homozygous variant position
						//
						contingencyTable[0][0] = coverage;
						contingencyTable[0][1] = 0;
						contingencyTable[1][0] = rawTumorAlleleCount[altLetter];
						contingencyTable[1][1] = rawTumorAlleleCount[refLetter];
						prob = kt_fisher_exact(contingencyTable[0][0], contingencyTable[0][1], contingencyTable[1][0], contingencyTable[1][1], &left, &right, &homTwotail);
						
						if(homTwotail>=GERMLINE_VARIANT_FET_CUTOFF) {
							// the test result is not significant, so I cannot reject the hypothesis that this position is homo-variant
							//
							tumorFormat = "1/1:";
						}
						else {
							tumorFormat = "0/1:";
						}
						tumorFormat += (IntToString(rawS[0].length()-rawTumorAlleleCount['?']) + ":" + IntToString(rawTumorAlleleCount[refLetter]) + "," + IntToString(rawTumorAlleleCount[altLetter]) + ":");
						
						if(rawTumorAlleleCount[refLetter]==0) {
							tumorFormat += "0,";
						}
						else {
							int tmp = round((double)rawTumorAlleleQuality[refLetter]/(double)rawTumorAlleleCount[refLetter]);
							tumorFormat += (IntToString(tmp) + ",");
						}
						if(rawTumorAlleleCount[altLetter]==0) {
							tumorFormat += "0:2";
						}
						else {
							int tmp = round((double)rawTumorAlleleQuality[altLetter]/(double)rawTumorAlleleCount[altLetter]);
							tumorFormat += (IntToString(tmp) + ":2");
						}
					}
				}
			}
			
			// output
			//
			outf << altString << '\t' << tumorFormat << '\t' << normalFormat << '\n';
		}
	}
	
	bam_mplp_destroy(iter);
	
	for(int i = 0; i < n; ++i) {
		if(data[i]->iter) {
			if(data[i]->iter->off) {
				free(data[i]->iter->off);
				data[i]->iter->off = NULL;
			}
			free(data[i]->iter);
			data[i]->iter = NULL;
		}
	}
	
	free(plp);
	free(n_plp);

	return 0;
}

void WriteHeader(std::ofstream &outf, int argc, char * const argv[], bam_header_t *tumorH, bam_header_t *normalH, const char *refGenome) {
	// muse version info
	//
	outf << "##MuSE_version=\"" << Version << " Build Date " << buildDate << " Build Time " << buildTime << "\"\n";
	// command line
	//
	outf << "##MuSE_call=\"";
	for(int i = 0; i < argc-1; i++) {
		outf << argv[i] << " ";
	}
	outf << argv[argc-1] << "\"\n";
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
			outf << "##TUMOR=\"Sample=" << kSM.s << ",File=" << argv[argc-2] << "\"\n";
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
			outf << "##NORMAL=\"Sample=" << kSM.s << ",File=" << argv[argc-1] << "\"\n";
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
			outf << "##contig=<ID=" << kSN.s << ",length=" << kLN.s;
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
				outf << ",assembly=" << kAS.s << ">\n";
				free(kAS.s);
			}
			else {
				outf << ">\n";
			}
		} else break;
		p = (sn>ln) ? ((sn>as)?sn:as) : ((ln>as)?ln:as);
	}
	// reference file
	//
	outf << "##reference=file://" << refGenome << '\n';
}

int MuseCall(int argc, char *argv[]) {
	// for dealing with program options, please refer to the way of PHYML
	//
	int         c;
    int         use_orphan = 0;
	const char *refGenome  = NULL;
    const char *regionFile = NULL;
    const char *outFile    = NULL;
	int         lineNumber = -1;
	mplp_conf_t mplp;
	memset(&mplp, 0, sizeof(mplp_conf_t));
	mplp.reg               = NULL;
	mplp.max_depth         = 8000;
	mplp.depth_cutoff      = 1;
	mplp.minAltFraction    = 0.005;
	mplp.min_output_brlens = 1e-4;
	mplp.flag              = MPLP_NO_ORPHAN | MPLP_REALN;
	mplp.flag             &= ~MPLP_REALN;
	while((c = getopt(argc, argv, "f:r:l:O:")) >= 0) {
		switch(c) {
			case 'f': refGenome           = optarg;         break;
			case 'r': mplp.reg            = strdup(optarg); break;
			case 'l': regionFile          = optarg;         break;
			case 'O': outFile             = optarg;         break;
		}
	}
	if(use_orphan)
		mplp.flag &= ~MPLP_NO_ORPHAN;
	if(argc == 1) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   MuSE call [options] tumor.bam matched_normal.bam\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "         -f FILE    faidx indexed reference sequence file\n");
		fprintf(stderr, "         -r STR     single region (chr:pos-pos) where somatic\n");
		fprintf(stderr, "                    mutations are called\n");
		fprintf(stderr, "         -l FILE    list of regions (chr:pos-pos or BED), one\n");
		fprintf(stderr, "                    region per line\n");
		fprintf(stderr, "         -O STR     output file name (suffix '.MuSE.txt' is\n");
		fprintf(stderr, "                    automatically added)\n");
		fprintf(stderr, "\n");
		return 1;
	}
	
	// check reference genome and output file
	//
	if(!refGenome) {
		std::cerr << "Please designate the reference genome fasta file.\n";
		return 0;
	}
	if(!outFile) {
		std::cerr << "Please designate the output file name.\n";
		return 0;
	}
	mplp.fai = fai_load(refGenome);
	if(mplp.fai == 0) {
		std::cerr << "Failed to load the reference genome index.\n";
		return 1;
	}

	// preparation ...
	//
	mplp_aux_t **data = (mplp_aux_t**)calloc(2, sizeof(void*)); // 0 for tumor; 1 for normal
	data[0]           = (mplp_aux_t*)calloc(1, sizeof(mplp_aux_t));
	data[1]           = (mplp_aux_t*)calloc(1, sizeof(mplp_aux_t));
	data[0]->fp       = bam_open((argv + optind)[0], "r");
	if(data[0]->fp == 0) {
		fprintf(stderr, "[%s] ERROR: fail to open file '%s'.\n", __func__, (argv + optind)[0]);
		return -1;
	}
	data[1]->fp = bam_open((argv + optind)[1], "r");
	if(data[1]->fp == 0) {
		fprintf(stderr, "[%s] ERROR: fail to open file '%s'.\n", __func__, (argv + optind)[1]);
		return -1;
	}
	data[0]->conf = &mplp;
	data[1]->conf = &mplp;
	data[0]->h    = bam_header_read(data[0]->fp);
	data[1]->h    = bam_header_read(data[1]->fp);
	
	bam_sample_t *sm = bam_smpl_init();
	bam_smpl_add(sm, (argv + optind)[0], (mplp.flag&MPLP_IGNORE_RG) ? 0 : data[0]->h->text);
	bam_smpl_add(sm, (argv + optind)[1], (mplp.flag&MPLP_IGNORE_RG) ? 0 : data[1]->h->text);

	if(mplp.max_depth * sm->n > 1<<20)
		fprintf(stderr, "(%s) Max depth is above 1M. Potential memory hog!\n", __func__);

	bam_index_t *idxTumor  = bam_index_load((argv + optind)[0]);
	bam_index_t *idxNormal = bam_index_load((argv + optind)[1]);
	if(idxTumor==0 || idxNormal==0) {
		fprintf(stderr, "[%s] fail to load index for the input.\n", __func__);
		exit(1);
	}

	char *ref = NULL;
	int  ref_len, ref_tid = -1;
	for(int i = 0; i < argc - optind; ++i)
		data[i]->ref_ptr = &ref, data[i]->ref_id = &ref_tid, data[i]->ref_len = &ref_len;
	
	// start ...
	//
	if(mplp.reg && regionFile) {
		fprintf(stderr, "[%s] option -l and -r are mutual exclusive.\n", __func__);
		exit(1);
	}

	std::string fileName(outFile);
	fileName += ".MuSE.txt";

	if(mplp.reg) {
		// open an output file
		//
		std::ofstream outf(fileName.c_str());
		WriteHeader(outf, argc, argv, data[0]->h, data[1]->h, refGenome);
		outf.precision(5);
		
		// pileup
		//
		mpileup(&mplp, argc - optind, data, idxTumor, idxNormal, &ref, ref_tid, ref_len, outf);
		outf.close();
		free(mplp.reg);
	}
	else if(regionFile) {
		std::ifstream regionf(regionFile);
		if(!regionf) {
			std::cerr << "Cannot open the region file.\n";
			exit(-1);
		}
		else {
			// open an output file
			//
			std::ofstream outf(fileName.c_str());
			WriteHeader(outf, argc, argv, data[0]->h, data[1]->h, refGenome);
			outf.precision(5);

			std::string word;
			regionf >> word;
			
			if(word.find(':') != std::string::npos && word.find('-') != std::string::npos) {
				std::vector<std::string> regions;
				regions.push_back(word);
				while(regionf >> word)
					regions.push_back(word);
				if(lineNumber >= 0 && lineNumber < regions.size()) {
					// pileup
					//
					mplp.reg = strdup(regions[lineNumber].c_str());
					mpileup(&mplp, argc - optind, data, idxTumor, idxNormal, &ref, ref_tid, ref_len, outf);
					free(mplp.reg);
				}
				else {
					for(unsigned r = 0; r < regions.size(); r++) {
						// pileup
						//
						mplp.reg = strdup(regions[r].c_str());
						mpileup(&mplp, argc - optind, data, idxTumor, idxNormal, &ref, ref_tid, ref_len, outf);
						free(mplp.reg);
					}
				}
			}
			else {
				std::vector<std::string> regions;
				regions.push_back(word);
				while(regionf >> word)
					regions.push_back(word);
				if(regions.size() % 3 != 0) {
					std::cerr << "The BED file has a wrong format.\n";
					exit(-1);
				}
				if(lineNumber >= 0 && lineNumber < regions.size()/3) {
					// pileup
					//
					mplp.reg = strdup((regions[lineNumber*3] + ":" + IntToString(StringToInt(regions[lineNumber*3+1])+1) + "-" + regions[lineNumber*3+2]).c_str());
					mpileup(&mplp, argc - optind, data, idxTumor, idxNormal, &ref, ref_tid, ref_len, outf);
					free(mplp.reg);
				}
				else {
					for(unsigned r = 0; r < regions.size(); r += 3) {
						// pileup
						//
						mplp.reg = strdup((regions[r] + ":" + IntToString(StringToInt(regions[r+1])+1) + "-" + regions[r+2]).c_str());
						mpileup(&mplp, argc - optind, data, idxTumor, idxNormal, &ref, ref_tid, ref_len, outf);
						free(mplp.reg);
					}
				}
			}
			outf.close();
			regionf.close();
		}
	}
	else {
		// without mplp.reg and regionFile, the default is analyzing all possible positions in the BAMs
		// open an output file
		//
		std::ofstream outf(fileName.c_str());
		WriteHeader(outf, argc, argv, data[0]->h, data[1]->h, refGenome);
		outf.precision(5);
		mpileup(&mplp, argc - optind, data, idxTumor, idxNormal, &ref, ref_tid, ref_len, outf);
		outf.close();
		if(mplp.reg != NULL)
			free(mplp.reg);
	}
	
	// clean ...
	//
	if(mplp.fai)
		fai_destroy(mplp.fai);
	
	for(int i = 0; i < argc - optind; ++i) {
		bam_close(data[i]->fp);
		if(data[i]->h)
			bam_header_destroy(data[i]->h);
		if(data[i]->iter)
			free(data[i]->iter);
		free(data[i]);
	}
	free(data);

	bam_smpl_destroy(sm);
	bam_index_destroy(idxTumor);
	bam_index_destroy(idxNormal);
	
	if(ref)
		free(ref);

	return 0;
}
//================================================================================================= MuSE sump
double Beta_Lk() {
	beta_lnlike = 0.0;
	for(int i = 0; i < altf.size(); i++) {
		beta_lnlike += dbeta(altf[i], betaShape[0], betaShape[1], 1);
	}
	int exponent = (int)FLOOR(log10(FABS(beta_lnlike)));
	if(sizeof(double) == 4) {
		min_diff_beta_lk = POW(10.,exponent - FLT_DIG + 1);
	}
	if(sizeof(double) == 8) {
		min_diff_beta_lk = POW(10.,exponent - DBL_DIG + 1);
	}
	return beta_lnlike;
}

double Negative_Misclassification_Prob() {
	 negative_misclassification_prob = -(pnorm(gmmCutoff, finalMu[0], finalSigma[0], 0, 0)*finalLambda[0] + pnorm(gmmCutoff, finalMu[1], finalSigma[1], 1, 0)*finalLambda[1]);
	int exponent = (int)FLOOR(log10(FABS(negative_misclassification_prob)));
	if(sizeof(double) == 4) {
		min_diff_misclassification_prob = POW(10.,exponent - FLT_DIG + 1);
	}
	if(sizeof(double) == 8) {
		min_diff_misclassification_prob = POW(10.,exponent - DBL_DIG + 1);
	}
	return negative_misclassification_prob;
}

/* function from mixtools v1.0.2
 * Compute the matrix of "posterior" probabilities in a finite mixture of univariate normal densities.
 * The algorithm used is fairly safe from a numerical perspective; it avoids over- or under-flow
 * as long as the values of sigma are not too large or small.
 */
void normpost(
			  int nn, /* sample size */
			  double *data,  /* vector of observations */
			  double *mu, /* current vector of component means */
			  double *sigma, /* current vector of component stdevs */
			  double *lambda, /* current vector of mixing parameters */
			  double *res2, /* n by m matrix of squared residuals */
			  double *work, /* 3*m-vector of workspace, which will be broken into 3 parts */
			  double *post, /* n by m matrix of posterior probabilities */
			  double *loglik /* scalar loglikelihood value */
) {
	int n=nn, m=2, i, j, minj=0;
	double x, r, rowsum, min=0.0;
	double *LamSigRatio = work+m; /* Second 1/3 of workspace, for frequently used constants */
	double *logLamSigRatio = work+2*m; /* Third 1/3 of workspace, for frequently used constants */
	
	*loglik = -(double)n * 0.91893853320467274178; /* n/2 times log(2pi) */
	for (j=0; j<m; j++){ /* store some often-used values to save time later */
		LamSigRatio[j] = lambda[j] / sigma[j];
		logLamSigRatio[j] = log(LamSigRatio[j]);
	}
	for (i=0; i<n; i++){
		x=data[i];
		for (j=0; j<m; j++) {
			r = x-mu[j];
			r = r*r;
			res2[i + n*j] = r;
			work[j] = r = r / (2.0 * sigma[j] * sigma[j]);
			/* Keep track of the smallest standardized squared residual.
			 By dividing everything by the component density with the
			 smallest such residual, the denominator of the posterior
			 is guaranteed to be at least one and cannot be infinite unless
			 the values of lambda or sigma are very large or small.  This helps
			 prevent numerical problems when calculating the posteriors.*/
			if (j==0 || r < min) {
				minj = j;
				min = r;
			}
		}
		/* At this stage, work contains the squared st'dized resids over 2 */
		rowsum = 1.0;
		for (j=0; j<m; j++) {
			if (j==minj)
				work[j] = 1.0;
			else {
				work[j] = (LamSigRatio[j] / LamSigRatio[minj]) * exp(min - work[j]);
				rowsum += work[j];
			}
		}
		/* At this stage, work contains the normal density at data[i]
		 divided by the normal density with the largest st'dized resid
		 Thus, dividing by rowsum gives the posteriors: */
		for (j=0; j<m; j++) {
			post[i + n*j] = work[j] / rowsum;
		}
		/* Finally, adjust the loglikelihood correctly */
		*loglik += log(rowsum) - min + logLamSigRatio[minj];
	}
}

int MuseSump(int argc, char *argv[]) {
	int        c;
	int        minDepth     = 8;
	double     passVAF      = 0.02;
	double     baseVAF      = 0.01;
	double     vafRatio     = 0.05;
	double     minTumorVAF  = 0.005;
	const char *outFile     = NULL;
	const char *inFile      = NULL;
	const char *dbsnpFile   = NULL;
	bool       isWGS        = false;
	bool       isWES        = false;
	
	// command options
	//
	while((c = getopt(argc, argv, "I:O:D:GE")) >= 0) {
		switch(c) {
			case 'I': inFile    = optarg; break;
			case 'O': outFile   = optarg; break;
			case 'D': dbsnpFile = optarg; break;
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
		fprintf(stderr, "         -D FILE    dbSNP vcf file that should be bgzip compressed,\n");
		fprintf(stderr, "                    tabix indexed and based on the same reference\n");
		fprintf(stderr, "                    genome used in 'MuSE call'\n");
		fprintf(stderr, "\n");
		return 1;
	}

	// must identify the data type but cannot be both
	//
	if(isWGS && isWES) {
		std::cerr << "Option -G and -E cannot be selected at the same time.\n";
		return 1;
	}
	if((!isWGS) && (!isWES)) {
		std::cerr << "Must identify the sequencing data type using either option -G or option -E.\n";
		return 1;
	}

	// check input and output files
	//
	if(!inFile) {
		std::cerr << "Please designate the input file.\n";
		return 0;
	}
	if(!outFile) {
		std::cerr << "Please designate the output file name.\n";
		return 0;
	}
	
	// check if dbSNP file was bgzipped
	//
	if(dbsnpFile && bgzf_is_bgzf(dbsnpFile)!=1) {
		fprintf(stderr,"Was bgzip used to compress %s?\n", dbsnpFile);
		return 1;
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
			return 1;
		}
		free(fnidx);
	}
	
	tabix_t *t;
	ti_iter_t iter;
	if(dbsnpFile) {
		if((t = ti_open(dbsnpFile, 0)) == 0) {
			fprintf(stderr, "Failed to open %s.\n", dbsnpFile);
			return 1;
		}
		if(ti_lazy_index_load(t) < 0) {
			fprintf(stderr,"Failed to load the index of %s.\n", dbsnpFile);
			return 1;
		}
	}
	
	// start ...
	//
	std::ifstream inf;
	inf.open(inFile);
	if(!inf) {
		std::cerr << "Failed to open " << inFile << '\n';
		return 0;
	}
	
	// load the input file content and remove unused lines
	//
	std::vector<std::string> Chrm;
	std::vector<std::string> Position;
	std::vector<std::string> Ref;
	std::vector<std::string> Alt;
	std::vector<double>      TaltF;
	std::vector<int>         ReadGroupCount;
	std::vector<bool>        IsDBSNP;
	std::vector<std::string> RSID;
	std::vector<std::string> Header;
	std::vector<bool>        UseSecondAllele;
	std::vector<int>         TumorRefCount;
	std::vector<int>         TumorAltCount;
	std::vector<int>         NormalRefCount;
	std::vector<int>         NormalAltCount;
	std::vector<std::string> AltString;
	std::vector<std::string> TumorFormat;
	std::vector<std::string> NormalFormat;
	
	Chrm.reserve(1000000);
	Position.reserve(1000000);
	Ref.reserve(1000000);
	Alt.reserve(1000000);
	TaltF.reserve(1000000);
	ReadGroupCount.reserve(1000000);
	IsDBSNP.reserve(1000000);
	RSID.reserve(1000000);
	Header.reserve(1000000);
	UseSecondAllele.reserve(1000000);
	TumorRefCount.reserve(1000000);
	TumorAltCount.reserve(1000000);
	NormalRefCount.reserve(1000000);
	NormalAltCount.reserve(1000000);
	AltString.reserve(1000000);
	TumorFormat.reserve(1000000);
	NormalFormat.reserve(1000000);
	
	std::string line = "";
	while(!inf.eof()) {
		getline(inf, line);
		if(line != "") {
			if(line[0]=='#') {
				Header.push_back(line);
				continue;
			}
			std::istringstream iss(line);
			std::vector<std::string> elements((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());
			if(elements.size() != INPUT_COLUMN_COUNT) {
				std::cerr << "Incomplete line of the input file.\n";
				inf.close();
				return 0;
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

						Chrm.push_back(elements[0]);
						Position.push_back(elements[1]);
						Ref.push_back(elements[2]);
						Alt.push_back(elements[3]);
						TaltF.push_back(tAltF);
						ReadGroupCount.push_back(nRG);
						TumorRefCount.push_back(tRefCount);
						TumorAltCount.push_back(tAltCount);
						NormalRefCount.push_back(nRefCount);
						NormalAltCount.push_back(nAltCount);
						if(elements[22]=="Y")
							UseSecondAllele.push_back(true);
						else
							UseSecondAllele.push_back(false);
						
						std::string rsID = "NA";
						if(dbsnpFile) {
							std::string coordinate = elements[0]+":"+elements[1]+"-"+elements[1];
							const char *s;
							int tid, beg, end, len;
							if(ti_parse_region(t->idx, coordinate.c_str(), &tid, &beg, &end) == 0) {
								iter = ti_queryi(t, tid, beg, end);
								if((s=ti_read(t, iter, &len)) != 0) {
									std::string dbsnpInfo(s);
									std::istringstream tmpISS(dbsnpInfo);
									std::vector<std::string> tmpElements((std::istream_iterator<std::string>(tmpISS)), std::istream_iterator<std::string>());
									rsID = tmpElements[2];
									IsDBSNP.push_back(true);
									RSID.push_back(rsID);
								}
								else {
									IsDBSNP.push_back(false);
									RSID.push_back(rsID);
								}
								ti_iter_destroy(iter);
							}
							else {
								IsDBSNP.push_back(false);
								RSID.push_back(rsID);
							}
						}
						else {
							IsDBSNP.push_back(false);
							RSID.push_back(rsID);
						}
						
						AltString.push_back(elements[27]);
						TumorFormat.push_back(elements[28]);
						NormalFormat.push_back(elements[29]);
					}
				}
			}
		}
	}
	inf.close();
	if(dbsnpFile)
		ti_close(t);

	// check if the BAM has multiple read group
	//
	bool isMultiReadGroup = false;
	for(int i = 0; i < ReadGroupCount.size(); i++) {
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
	for(int i = 0; i < Header.size(); i++) {
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
		return 1;
	}
	if(TumorSample.size()!=1) {
		std::cerr << "Either tumor sample ID is missing or multiple tumor sample IDs are found. Please check input file.\n";
		return 1;
	}
	if(NormalSample.size()!=1) {
		std::cerr << "Either normal sample ID is missing or multiple normal sample IDs are found. Please check input file.\n";
		return 1;
	}
	if(ReferenceGenome.size()>1) {
		std::cerr << "Multiple reference genomes used. Please check input file.\n";
		return 1;
	}
	
	// extract contig names, used for sorting "MuSE call" output
	//
	std::vector<std::string> ContigName;
	for(int i = 0; i < ContigInfo.size(); i++) {
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
	
	FinalChrm.reserve(1000000);
	FinalPosition.reserve(1000000);
	FinalRef.reserve(1000000);
	FinalAlt.reserve(1000000);
	FinalTaltF.reserve(1000000);
	FinalReadGroupCount.reserve(1000000);
	FinalIsDBSNP.reserve(1000000);
	FinalRSID.reserve(1000000);
	FinalUseSecondAllele.reserve(1000000);
	FinalTumorRefCount.reserve(1000000);
	FinalTumorAltCount.reserve(1000000);
	FinalNormalRefCount.reserve(1000000);
	FinalNormalAltCount.reserve(1000000);
	FinalAltString.reserve(1000000);
	FinalTumorFormat.reserve(1000000);
	FinalNormalFormat.reserve(1000000);

	typedef std::map<int, callInfo>::value_type callType;

	for(int i = 0; i < ContigName.size(); i++) {
		std::map<int, callInfo> ContigCall;
		std::vector<int> TmpPosition;
		TmpPosition.reserve(100000);
		for(int j = 0; j < Chrm.size(); j++) {
			if(Chrm[j]==ContigName[i]) {
				int position = atoi(Position[j].c_str());
				TmpPosition.push_back(position);
				
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
				
				ContigCall.insert(callType(position, oneCall));
			}
		}
		std::sort(TmpPosition.begin(), TmpPosition.end());
		for(int j = 0; j < TmpPosition.size(); j++) {
			FinalChrm.push_back(ContigName[i]);
			FinalPosition.push_back(TmpPosition[j]);
			FinalRef.push_back(ContigCall[TmpPosition[j]].ref);
			FinalAlt.push_back(ContigCall[TmpPosition[j]].alt);
			FinalTaltF.push_back(ContigCall[TmpPosition[j]].tAltF);
			FinalReadGroupCount.push_back(ContigCall[TmpPosition[j]].readGroupCount);
			FinalIsDBSNP.push_back(ContigCall[TmpPosition[j]].isDBSNP);
			FinalRSID.push_back(ContigCall[TmpPosition[j]].rsID);
			FinalUseSecondAllele.push_back(ContigCall[TmpPosition[j]].useSecondAllele);
			FinalTumorRefCount.push_back(ContigCall[TmpPosition[j]].tumorRefCount);
			FinalTumorAltCount.push_back(ContigCall[TmpPosition[j]].tumorAltCount);
			FinalNormalRefCount.push_back(ContigCall[TmpPosition[j]].normalRefCount);
			FinalNormalAltCount.push_back(ContigCall[TmpPosition[j]].normalAltCount);
			FinalAltString.push_back(ContigCall[TmpPosition[j]].altString);
			FinalTumorFormat.push_back(ContigCall[TmpPosition[j]].tumorFormat);
			FinalNormalFormat.push_back(ContigCall[TmpPosition[j]].normalFormat);
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
		lnTEVAF.reserve(FinalTaltF.size());
		if(isMultiReadGroup) {
			for(int i = 0; i < FinalTaltF.size(); i++) {
				if(FinalReadGroupCount[i] > 1) {
					lnTEVAF.push_back(log(FinalTaltF[i]));
				}
			}
		}
		else {
			for(int i = 0; i < FinalTaltF.size(); i++) {
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
		for(int i = 0; i <= maxIndexLeftPartition; i++) {
			leftPartitionMean += lnTEVAF[i];
		}
		leftPartitionMean = leftPartitionMean/(double)(maxIndexLeftPartition+1);
		double rightPartitionMean = 0.0;
		for(int i = minIndexRightPartition; i < dataSize; i++) {
			rightPartitionMean += lnTEVAF[i];
		}
		rightPartitionMean = rightPartitionMean/(double)(dataSize-minIndexRightPartition);
		double lnTEVAFMean = 0.0;
		for(int i = 0; i < dataSize; i++) {
			lnTEVAFMean += lnTEVAF[i];
		}
		lnTEVAFMean = lnTEVAFMean/(double)dataSize;
		// standard deviation
		//
		double leftPartitionSD = 0.0;
		for(int i = 0; i <= maxIndexLeftPartition; i++) {
			leftPartitionSD += (lnTEVAF[i] - leftPartitionMean)*(lnTEVAF[i] - leftPartitionMean);
		}
		leftPartitionSD = sqrt(leftPartitionSD/(double)(maxIndexLeftPartition));
		double originalLeftPartitionSD = leftPartitionSD;
		double rightPartitionSD = 0.0;
		for(int i = minIndexRightPartition; i < dataSize; i++) {
			rightPartitionSD += (lnTEVAF[i] - rightPartitionMean)*(lnTEVAF[i] - rightPartitionMean);
		}
		rightPartitionSD = sqrt(rightPartitionSD/(double)(dataSize-minIndexRightPartition-1));
		double originalRightPartitionSD = rightPartitionSD;
		double lnTEVAFSD = 0.0;
		for(int i = 0; i < dataSize; i++) {
			lnTEVAFSD += (lnTEVAF[i] - lnTEVAFMean)*(lnTEVAF[i] - lnTEVAFMean);
		}
		lnTEVAFSD = sqrt(lnTEVAFSD/(double)(dataSize-1));
		
		// try EM 50 times
		//
		EMEstimate.clear();
		for(int ithReplicate = 0; ithReplicate < emReplicate; ithReplicate++) {
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
					for(int i = 0; i < 2; i++) {
						lambda[i] = 0.0;
						for(int j = 0; j < dataSize; j++) {
							lambda[i] += postProbs[i*dataSize+j];
						}
						lambda[i] = lambda[i]/(double)dataSize;
					}
					for(int i = 0; i < 2; i++) {
						mu[i] = 0.0;
						for(int j = 0; j < dataSize; j++) {
							mu[i] += postProbs[i*dataSize+j] * lnTEVAF[j];
						}
						mu[i] = mu[i]/((double)dataSize*lambda[i]);
					}
					for(int i = 0; i < 2; i++) {
						sigma[i] = 0.0;
						for(int j = 0; j < dataSize; j++) {
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
							exit(-1);
						}
						break;
					}
					// E-step
					//
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
			if(mu[0]==mu[1])
				continue;
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
			//
			if(mu[1]>-4.605170185988091) {
				emEst tmpEstimate;
				tmpEstimate.loglik  = loglik;
				tmpEstimate.mu1     = mu[0];
				tmpEstimate.mu2     = mu[1];
				tmpEstimate.sigma1  = sigma[0];
				tmpEstimate.sigma2  = sigma[1];
				tmpEstimate.lambda1 = lambda[0];
				tmpEstimate.lambda2 = lambda[1];
				EMEstimate.push_back(tmpEstimate);
			}
		}
		
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
			for(int i = 0; i < EMEstimate.size(); i++) {
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
					for(int i = 0; i < 2; i++) {
						finalLambda[i] = 0.0;
						for(int j = 0; j < dataSize; j++) {
							finalLambda[i] += postProbs[i*dataSize+j];
						}
						finalLambda[i] = finalLambda[i]/(double)dataSize;
					}
					for(int i = 0; i < 2; i++) {
						finalSigma[i] = 0.0;
						for(int j = 0; j < dataSize; j++) {
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
							exit(-1);
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
				if(FABS(prob_new - prob_old) < min_diff_misclassification_prob)
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
		altf.reserve(FinalTaltF.size());
		for(int i = 0; i < FinalTaltF.size(); i++) {
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
			int i;
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
			if(FABS(lk_new - lk_old) < min_diff_beta_lk)
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
		for(int i = 0; i < EMEstimate.size(); i++) {
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
	for(int i = 0; i < MuSEVersion.size(); i++) {
		outf << MuSEVersion[i] << '\n';
	}
	for(int i = 0; i < MuSECallCMD.size(); i++) {
		std::string base = MuSECallCMD[i];
		std::string substr = "##MuSE_call_" + IntToString(i+1);
		base.replace(0, 11, substr);
		outf << base << '\n';
	}
	// sump command line
	//
	outf << "##MuSE_sump=\"";
	for(int i = 0; i < argc-1; i++) {
		outf << argv[i] << " ";
	}
	outf << argv[argc-1] << "\"\n";
	for(int i = 0; i < ContigInfo.size(); i++) {
		outf << ContigInfo[i] << '\n';
	}
	for(int i = 0; i < ReferenceGenome.size(); i++) {
		outf << ReferenceGenome[i] << '\n';
	}
	outf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL\n";
	std::string whichTier = "";
	for(int i = 0; i < FinalTaltF.size(); i++) {
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
	
	return 0;
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
	fprintf(stderr, "                 Raw output file will be created. If WGS\n");
	fprintf(stderr, "                 data is used, it is recommended to split\n");
	fprintf(stderr, "                 the WGS data into small blocks (<50Mb).\n");
	fprintf(stderr, "         sump    Generate final calls in VCF format.\n");
	fprintf(stderr, "                 Besides PASS calls, tier-based calls are\n");
	fprintf(stderr, "                 reported using dynamic cutoff searching\n");
	fprintf(stderr, "                 approach. If there are multiple raw output\n");
	fprintf(stderr, "                 files from WGS data, please have them\n");
	fprintf(stderr, "                 concatenated first.\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[]) {
	if(argc < 2)
		return usage();
	if(strcmp(argv[1], "call") == 0)      return MuseCall(argc-1, argv+1);
	else if(strcmp(argv[1], "sump") == 0) return MuseSump(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	return 0;
}







