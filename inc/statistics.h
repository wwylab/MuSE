#include <iostream>
#include <fstream>
#include <cstring>
#include <sys/stat.h>
#include <math.h>
#include <float.h>
#include "omp.h"
#include "muse_const.h"

using namespace std;


//================================================================================================= Beta Distribution
#define MY_PI          3.141592653589793238462643383280	/* pi */
#define M_2PI          6.283185307179586476925286766559	/* 2*pi */
#define M_LN_sqrt_2PI  0.918938533204672741780329736406	/* log(sqrt(2*pi)) == log(2*pi)/2 */
#define M_LN_sqrt_PId2 0.225791352644727432363097614947	/* log(sqrt(pi/2)) */
#define M_LOG10_2      0.301029995663981195213738894724	/* log10(2) */
#define MY_LN2         0.693147180559945309417232121458	/* ln(2) */
#define M_1_sqrt_2PI   0.398942280401432677939946059934	/* 1/sqrt(2pi) */
#define M_sqrt_32      5.656854249492380195206754896838	/* sqrt(32) */

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
#define R_D_LExp_toms708(x)   (log_p ? R_Log1_Exp_toms708(x) : r_log1p(-x))					/* MM added R_D_LExp, so redefine here in terms of rr_expm1 */
#define R_Log1_Exp_toms708(x) ((x) > -MY_LN2 ? log(-rr_expm1(x)) : r_log1p(-exp(x)))
#define R_D_Lval(p)				(lower_tail ? (p)				: (0.5 - (p) + 0.5)	)			/* Use 0.5 - p + 0.5 to perhaps gain 1 bit of accuracy: p     */
#define R_D_Cval(p)				(lower_tail ? (0.5 - (p) + 0.5) : (p)				)			/* Use 0.5 - p + 0.5 to perhaps gain 1 bit of accuracy: 1 - p */
#define R_DT_CIv(p)           (log_p ? (lower_tail ? -r_expm1(p): exp(p))     : R_D_Cval(p))	/* #define R_DT_CIv(p)	R_D_Cval(R_D_qIv(p))	*  1 - p in qF	*/

//#define ME_NONE      0	/*	no error */
#define ME_DOMAIN    1	/*	argument out of domain */
#define ME_RANGE     2	/*	value out of range */
#define ME_NOCONV    4	/*	process did not converge */
#define ME_PRECISION 8	/*	does not have "full" precision */
#define ME_UNDERFLOW 16	/*	and underflow occured (important for IEEE)*/

#define MATHLIB_ERROR(fmt,x)   { printf(fmt, x); exit(EXIT_FAILURE); }
#define MATHLIB_WARNING(fmt,x) printf(fmt,x)
#define ML_ERR_return_NAN      { ML_ERROR(ME_DOMAIN, ""); return ML_NAN; }
#define ML_ERROR(x, s) {										\
if(x > ME_DOMAIN) {												\
	string msg = "";												\
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
	MATHLIB_WARNING(msg.c_str(), s);									\
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

// for beta fit
//
extern std::vector<double> altf;
extern std::vector<double> betaShape;
extern double min_diff_beta_lk;
extern double beta_lnlike;

// for GMM fit
//
extern double finalMu[2];
extern double finalSigma[2];
extern double finalLambda[2];
extern double gmmCutoff;
extern double min_diff_misclassification_prob;
extern double negative_misclassification_prob;


// Define statistics functions
double unif_rand(void);

double exp_rand(void);

double dnorm4(double x, double mu, double sigma, int give_log);

void pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p);

double pnorm5(double x, double mu, double sigma, int lower_tail, int log_p);

double qnorm5(double p, double mu, double sigma, int lower_tail, int log_p);

double norm_rand(void);


double r_log1p(double x);


static double rlog1(double x);

double r_expm1(double x);

double dbeta(double x, double a, double b, int give_log);

double qbeta(double alpha, double p, double q, int lower_tail, int log_p);

double Generic_Brent_Lk(double *param, double ax, double cx, double tol, int n_iter_max, int quickdirty, double(*obj_func)());

void normpost(int nn, double *data, double *mu, double *sigma, double *lambda, double *res2, double *work, double *post, double *loglik);

double Beta_Lk();

double Negative_Misclassification_Prob();
