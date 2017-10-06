#include	<stdio.h>
#include	<stdlib.h>

#define         DEBUG(x)   printf("DEBUG %d\n", x); fflush(stdout);
#define         DEBUGi(x)   printf("DEBUG %d\n", x); fflush(stdout);
#define         DEBUGf(x)   printf("DEBUG %f\n", x); fflush(stdout);
#define         DEBUGs(x)   printf("DEBUG %s\n", x); fflush(stdout);
#define         SHOW_DIRS   printf("\n\nHERE IS YOUR 'DIR' ENVIRONMENT\n\n"); system("env | grep DIR\n");

#define	PRAHA_SAMPS_PER_SEC	10
#define	DEFAULT_DISTANCE_dT	5
#define	MAX_XSZ			21
#define	MAX_YSZ			21
#define	REDUCE_X		12.0
#define	REDUCE_Y		12.0
#define	MAX_UNITS		10
#define	MAX_FRAMES		2
#define	N_TRACK_HEADER_LINES	1
#define	MAX_COLORS		9
#define	WHITE			0
#define	YELLOW			1
#define	ORANGE			2
#define	RED			3
#define	GREEN			4
#define	BLUE			5
#define	PURPLE			6
#define	GRAY			7
#define	BLACK			8

#define         TWO_PI          6.2831853
#define         DEG_PER_RAD     57.29578
#define         RAD_PER_DEG     0.017453
#define		SMALL_NUMBER	0.00001

#define		N_THRESHOLD_TYPES	1
#define		THRESHOLD_CONSTANT	0.5
#define		MEAN		0
#define		CONSTANT	1
#define		ZERO		2
#define		STAT		3
#define		N_SD_THRESHOLD	1.00

#define		MIN_FLD_SZ	4
#define		MIN_FLD_Y	2
#define		MIN_FLD_X	2
#define		MIN_SAMPS_4_GOOD	3

typedef	char	str256[256];
typedef	char	BYTE;

typedef int	INT_MAP[MAX_YSZ][MAX_XSZ];
typedef	float	FLOAT_MAP[MAX_YSZ][MAX_XSZ];

typedef struct {
        char    line1[256];
        char    line2[256];
        char    experimenter[256];
        char    cell_name[256];
        char    date[256];
        char    time[256];
        char    condition[256];
        int     n_maps;
        int     n_units;
        int     Ysz;
        int     Xsz;
        } HEADER;

typedef struct {
	float	max;
	float	medn;
	} KEY;

typedef struct {
	int		num;
	int		size;
	float		time;
	int		spks;
	float		rate;
	int		oof_act_pix;
	float		oof_time;
	int		oof_spks;
	float		oof_rate;
	float		cent_y;
	float		cent_x;
	float		cent_ang;
	float		cent_rad;
	float		cent_rate;
	float		cg_y;
	float		cg_x;
	float		cg_ang;
	float		cg_rad;
	float		cg_rate;
	float		coh;
	float		concent;
	float		s2n;
	} FIELD;

typedef	struct {
	int		num;
	int		map;
	int		b_spks;
	int		brsts;
	int		s_spks;
	int		n_colors;
	int		area;
	float		tot_time;
	int		tot_spks;
	int		top_y;
	int		bot_y;
	int		rgt_x;
	int		lft_x;
	float		cent_x;
	float		cent_y;
	float		grand_rate;
	float		if_rate;
	float		if_time;
	int		if_spks;
	float		oof_rate;
	float		oof_time;
	int		oof_spks;
	float		in2out;
	float		rate_time_corr;
	int		n_rate_time_corr;
	int		act_pix;
	int		patches;
	float		z_coher;
	int		n_in_coh;
	float		info_cnt;
	float		z_info;
	float		threshold;
	int		n_fields;
	FIELD		*field;
	INT_MAP		time_map;
	INT_MAP		spk_map;
	INT_MAP		color_map;
	KEY		key[MAX_COLORS];
	FLOAT_MAP	rate_map;
	INT_MAP		field_map;
	float		samps_per_sec;
	int		*isi;
	int		isi_n;
	int		*psth;
	int		psth_n;
	int		*a_psdh;
	int		*d_psdh;
	char		file[256];
	FILE		*afp;
	FILE		*dfp;
	} FRAME_MAP;

typedef struct {
	FRAME_MAP	*fm;
	} UNIT;

typedef struct {
        int     f_nr;
        int     f_sz;
	} SORT_FLDS;

typedef struct {
	int	x;
	int	y;
	} POSITION;

typedef  struct {
	POSITION *p;
	} POSITION_LIST;

typedef int SPK_CNT;

typedef  struct{
	int		*s;
	unsigned long 	*ts;
	int		n;
	} UNIT_LIST;

typedef struct {
	POSITION_LIST	rf[MAX_FRAMES];
	UNIT_LIST	unit[MAX_UNITS];
	} TIME_SERIES;
