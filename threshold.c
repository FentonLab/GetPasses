#include <math.h>
#include "tsf.h"

float	get_stat_threshold(r, t, mean, n_sd, Ysz, Xsz)
	FLOAT_MAP	r;
	INT_MAP		t;
	float	mean, n_sd;
	int	Ysz, Xsz;
{
	int	y,x,n = 0;
	float	X, sum_XX = 0.0, sd, thresh;
	double	sqrt();

	for(y=0;y < Ysz;y++){
		for(x=0;x < Xsz;x++){
			if(t[y][x] > 0){
				X = r[y][x] - mean;
				sum_XX += (X * X);
				n++;
			}
		}
	}

	sd = sum_XX / (float)(n - 1);
        sd = sqrt(sd);

	thresh = mean + (n_sd * sd); 

	return(thresh);
}
