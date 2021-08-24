#ifndef __MUSE_CONST_H__
#define __MUSE_CONST_H__

#include "Reference.h"
#include <string>

/*! @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
#define BAM_FPAIRED        1
/*! @abstract the read is mapped in a proper pair */
#define BAM_FPROPER_PAIR   2
/*! @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
#define BAM_FUNMAP         4
/*! @abstract the mate is unmapped */
#define BAM_FMUNMAP        8
/*! @abstract the read is mapped to the reverse strand */
#define BAM_FREVERSE      16
/*! @abstract the mate is mapped to the reverse strand */
#define BAM_FMREVERSE     32
/*! @abstract this is read1 */
#define BAM_FREAD1        64
/*! @abstract this is read2 */
#define BAM_FREAD2       128
/*! @abstract not primary alignment */
#define BAM_FSECONDARY   256
/*! @abstract QC failure */
#define BAM_FQCFAIL      512
/*! @abstract optical or PCR duplicate */
#define BAM_FDUP        1024
/*! @abstract supplementary alignment */
#define BAM_FSUPP       2048


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
#define SIGN(a,b)      ((b) > 0.0 ? fabs(a) : -fabs(a))
#define For(i,n)       for(i=0; i<n; i++)
#define SHFT(a,b,c,d)  (a)=(b);(b)=(c);(c)=(d);
#define MAX(a,b)       ((a)>(b)?(a):(b))


#define BAM_DEF_MASK (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPP)

#define _cop(c) ((c)&BAM_CIGAR_MASK)
#define _cln(c) ((c)>>BAM_CIGAR_SHIFT)


class mplp_conf_t{
public:
	int     flag, flag_mask, max_depth, depth_cutoff;
	double  minAltFraction, min_output_brlens;
	Reference* ref;
	int argc;
	char** argv;
	mplp_conf_t(std::string& refName, int argc_in, char* argv_in[]): argc(argc_in), argv(argv_in){
		max_depth         = 8000;
		depth_cutoff      = 1;
		minAltFraction    = 0.005;
		min_output_brlens = 1e-4;
		flag              = MPLP_NO_ORPHAN;
		flag_mask = BAM_DEF_MASK;

		ref = new Reference();
		ref->openFile(refName.c_str());
		ref->readFile(false);
	}

	~mplp_conf_t(){
		delete ref;
	}
};

struct cstate_t{
	int k, x, y;
};

struct bam_pileup1_t_pb{
	bam1_t *b;
	int64_t end;
	int32_t qpos;
	int indel;
	cstate_t s;
	bool is_del, is_head, is_tail, is_refskip;
};



#endif