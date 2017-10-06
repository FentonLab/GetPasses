#include	<math.h>

float	correl_2dim_arr(arr,n)

float	arr[][2];
int	n;

{
	float	a,aa,b,bb,ab;
	float	t1,t2;
	double	sqrt();
	int	i;

	a = aa = b = bb = ab = 0.0;

	for(i=0;i<n;i++){
		t1 = arr[i][0];
		t2 = arr[i][1];
		a += t1;
		aa += t1*t1;
		b += t2;
		bb += t2*t2;
		ab += t1*t2;
		}

	t1 = ab - a*b/n;
	t2 = (aa - a*a/n)*(bb - b*b/n);
	return(t1/sqrt(t2));
}	

float	z_trans(r)
float	r;

{
	double	log();

	return(0.5*log((1.+ r)/(1. - r)));
}

float	inv_z_trans(z)
float	z;

{
	float	exp2z;
	double	exp();

	exp2z = exp(2.*z);

	return((exp2z - 1)/(exp2z+1));
}

float	t_of_correl(cor,n)
float	cor;
int	n;

{
	float	t1,t2;
	double	sqrt();

	t1 = cor*sqrt((double) n - 2.0);
	t2 = sqrt(1.0 - cor*cor);

	return(t1/t2);
}
