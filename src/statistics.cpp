#include "statistics.h"

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

typedef enum {
	BUGGY_KINDERMAN_RAMAGE, AHRENS_DIETER, BOX_MULLER, USER_NORM, INVERSION, KINDERMAN_RAMAGE
} N01type;


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
	return (give_log ? -(M_LN_sqrt_2PI + 0.5*x*x + log(sigma)) : M_1_sqrt_2PI * exp(-0.5*x*x) / sigma);
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
	else if (y <= M_sqrt_32) {
		
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
			 *cum = -.5*xsq - M_LN_sqrt_2PI - log(x) + r_log1p(-del);
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
		temp = (M_1_sqrt_2PI - temp) / y;
		
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
	/* R has  M_1_sqrt_2PI */
	
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
	/* R has  M_1_sqrt_2PI , and M_LN_sqrt_2PI = ln(sqrt(2*pi)) = 0.918938.. */
	
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
	
	return(log_p ? -M_LN_sqrt_2PI + .5*log(b * x0) + z - bcorr(a,b) : const__ * sqrt(b * x0) * z * exp(-bcorr(a, b)));
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
			return(M_LN_sqrt_2PI + (x - 0.5) * log(x) - x);
		else
			return M_LN_sqrt_2PI + (x - 0.5) * log(x) - x + lgammacor(x);
	}
	/* else: x < -10; y = -x */
	sinpiy = fabs(sin(MY_PI * y));
	
	if (sinpiy == 0) { /* Negative integer argument === Now UNNECESSARY: caught above */
		MATHLIB_WARNING(" ** should NEVER happen! *** [lgamma.c: Neg.int, y=%g]\n",y);
		ML_ERR_return_NAN;
	}
	
	ans = M_LN_sqrt_PId2 + (x - 0.5) * log(y) - x - log(sinpiy) - lgammacor(y);
	
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
		return(lgammafn(n + 1.) - (n + 0.5)*log(n) + n - M_LN_sqrt_2PI);
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
			value = exp((y - 0.5) * log(y) - y + M_LN_sqrt_2PI + ((2*y == (int)2*y)? stirlerr(y) : lgammacor(y)));
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
		return log(q) * -0.5 + M_LN_sqrt_2PI + corr + (p - 0.5) * log(p / (p + q)) + q * r_log1p(-p / (p + q));
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
		
		if((fabs(fu-old_lnL) < tol) || (iter > n_iter_max - 1)) {
			(*param) = x;
			fu = (*obj_func)();
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

//================================================================================================= MuSE sump

double Beta_Lk() {
    beta_lnlike = 0.0;
    for(int i = 0; i < altf.size(); i++) {
        beta_lnlike += dbeta(altf[i], betaShape[0], betaShape[1], 1);
    }
    int exponent = (int)floor(log10(fabs(beta_lnlike)));
    if(sizeof(double) == 4) {
        min_diff_beta_lk = pow(10.,exponent - FLT_DIG + 1);
    }
    if(sizeof(double) == 8) {
        min_diff_beta_lk = pow(10.,exponent - DBL_DIG + 1);
    }
    return beta_lnlike;
}

double Negative_Misclassification_Prob() {
	 negative_misclassification_prob = -(pnorm(gmmCutoff, finalMu[0], finalSigma[0], 0, 0)*finalLambda[0] + pnorm(gmmCutoff, finalMu[1], finalSigma[1], 1, 0)*finalLambda[1]);
	int exponent = (int)floor(log10(fabs(negative_misclassification_prob)));
	if(sizeof(double) == 4) {
		min_diff_misclassification_prob = pow(10.,exponent - FLT_DIG + 1);
	}
	if(sizeof(double) == 8) {
		min_diff_misclassification_prob = pow(10.,exponent - DBL_DIG + 1);
	}
	return negative_misclassification_prob;
}

