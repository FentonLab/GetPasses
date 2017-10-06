#include <math.h>
#include "tsf.h"

float	rt_correl(arr1, arr2, nin, nout)
float	arr1[], arr2[];
int	nin, *nout;
{
	float	a, aa, b, bb, ab;
	float	t1, t2;
	double	sqrt();
	int	i;

	*nout = a = aa = b = bb = ab = 0.0;

	for (i = 0; i < nin; i++) {
		t1 = arr1[i];
		t2 = arr2[i];
		if (t1 == -1.0 || t2 == -1.0)
			continue;
		(*nout)++;
		a += t1;
		aa += t1 * t1;
		b += t2;
		bb += t2 * t2;
		ab += t1 * t2;
	}

	t1 = ab - a * b / (*nout);
	t2 = (aa - a * a / (*nout)) * (bb - b * b / (*nout));
	if (t2 <= 0.0) {
		/* printf("%ccorrelation error: num = %f, denom = %f\n", 7, t1 , t2); */
		return(0.0);
	}

	return(t1 / sqrt(t2));
}


get_time_rate_corr(fm, t, r, Ysz, Xsz)
FRAME_MAP	*fm;
INT_MAP	t;
FLOAT_MAP r;
int	Ysz, Xsz;
{
	int	i, j;
	float	f[MAX_YSZ * MAX_XSZ], g[MAX_YSZ * MAX_XSZ];
	float	rt_correl();

	for (i = 0; i < Ysz; i++)
		for (j = 0; j < Xsz; j++) {
			f[i*Ysz + j] = (float) t[i][j];
			g[i*Ysz + j] = r[i][j];
		}

	fm->rate_time_corr = rt_correl(f, g, Ysz * Xsz, &(fm->n_rate_time_corr));
}
