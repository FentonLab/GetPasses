#define	SAMPS_PER_SEC		60
#define	SEC_PER_APPROACH	5
#define	SEC_PER_DEPARTURE	5
#define	N_APPROACH_SAMPS	(SAMPS_PER_SEC * SEC_PER_APPROACH)
#define	N_DEPARTURE_SAMPS	(SAMPS_PER_SEC * SEC_PER_DEPARTURE)
#define	DWELL_CRITERION		5
#define	TRACKER_RESOLUTION_cmpp 0.33  /* 0.7 */
#define	TARGET_RADIUS_cm	10.0 * 1.1
#define	TARGET_RADIUS_pix	(TARGET_RADIUS_cm / TRACKER_RESOLUTION_cmpp)

#define	APPARATUS_RADIUS_cm	40.0
#define	APPARATUS_RADIUS_pix	(APPARATUS_RADIUS_cm / TRACKER_RESOLUTION_cmpp)
#define CENTER_X		10.0
#define CENTER_Y		10.0
#define MAX_SPIKES		6000
#define MAX_ISI			2000
#define MAX_PSDH_BINS		255
#define PSDH_BIN_SZ		5	/* gives a pixel resolution i.e. 5 cm per bin */

#define	APPROACH		0
#define NEUTRAL			-1
#define	DEPARTURE		1

#define		HOT_SPOT	1000
#define		SAMPS_PER_FRAME	2
#define		MAX_SHORT      65535

typedef unsigned short int USI;

typedef struct{
	int	start;
	int	stop;
}EPISODE;

typedef struct {
	double	x;
	double	y;
	double	cx;
	double	cy;
	double	r;
	INT_MAP map;
	int	enabled;
} TARGET;
