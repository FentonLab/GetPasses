#include	<math.h>
#include	<sys/time.h>
#include	<stdio.h>
#include	"tsf.h"
#include	"my_math.h"

#define	INFO_ITERATIONS	100
#define MAX_SHIFT	(PRAHA_SAMPS_PER_SEC * 500)

float	get_info_cnt(rate,t,R,T,Ysz,Xsz)
FLOAT_MAP	rate;
INT_MAP	t;
float	R, T;
int	Ysz,Xsz;
{
	int	i,j;
	float	p,r;
	float ic = 0.0;

	for(i=0;i<Ysz;i++)
		for(j=0;j<Xsz;j++){
			if(t[i][j] <= 0)
				continue;
			p = (float)t[i][j] / T;
			r = rate[i][j] / R;
			if(r > 0)
				ic += r * (float)log2((double)r) * p;
		}
	return(ic);
}		

float	z_info_cnt(p, spks, n_samps, IC, Ysz, Xsz, samps_per_sec)
POSITION	*p;
int		*spks;
float		IC, samps_per_sec;
int		n_samps, Ysz, Xsz;
{
	float	R, T;
	double	drand48(), sqrt();
	void	srand48();
	int	gettimeofday();
	struct timeval tp;
	struct timezone tzp;

	static INT_MAP		t;
	static FLOAT_MAP	r;
	static INT_MAP		s;
	float	z, XX, sum, mean, sd, ic;
	int	x, y, i, j, k, shift, S, n = INFO_ITERATIONS + 1; 

	srand48((long)gettimeofday(&tp, &tzp));

	T = (float)n_samps;	
	sum = XX = 0.0;

	/* make a time map bc the tsf map doesn't exactly correspond to the track file */ 
	for(y = 0; y < Ysz; y++)
		for(x = 0; x < Xsz; x++)
			t[y][x] = 0;

	for(i = 0; i < n_samps; i++)
		t[p[i].y][p[i].x]++;

	for(i = 0; i < n; i++){
		shift = (int)(MAX_SHIFT * drand48());
		S = 0;
		for(y = 0; y < Ysz; y++){
			for(x = 0; x < Xsz; x++){
				s[y][x] = 0;
				r[y][x] = -1.0;
			}
		}
		for(j = 0; j < n_samps; j++){
			k = ((j + shift) < n_samps) ? (j + shift) : (j + shift - n_samps);
			s[p[j].y][p[j].x] += spks[k];
		}	

		for(y = 0; y < Ysz; y++){
			for(x = 0; x < Xsz; x++){
				if(t[y][x] > 0)
					r[y][x] = (float)s[y][x] / (float)t[y][x] * samps_per_sec;
					S += s[y][x];
			}
		}
		R = (float)(S / T * samps_per_sec);
		ic = get_info_cnt(r,t,R,T,Ysz,Xsz);
/*
printf("%d\t%f\t%f\t%f\t%d\t%f\t%d\n",shift, IC, ic, R, S, T, XXX);
*/
		sum += ic;
		XX += (ic * ic);
	}
/*
	for(y = 0; y < Ysz; y++){
		for(x = 0; x < Xsz; x++)
			printf("%d  ",t[y][x]);
		printf("\n");
	}
	for(y = 0; y < Ysz; y++){
		for(x = 0; x < Xsz; x++)
			printf("%d  ",s[y][x]);
		printf("\n");
	}
	for(y = 0; y < Ysz; y++){
		for(x = 0; x < Xsz; x++)
			printf("%3.1f  ",r[y][x]);
		printf("\n");
	}
*/

	mean = sum / (float)n;
	sd = (((float)n * XX) - (sum * sum)) / (n * (n - 1));
	sd = (float)sqrt((double)sd);
	
	z = (IC - mean) / sd;
		
	return (z);
}
