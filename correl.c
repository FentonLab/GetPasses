#include <math.h>

float	correl(arr1,arr2,nin,nout)

float	arr1[],arr2[];
int	nin,*nout;

{
	float	a,aa,b,bb,ab;
	float	t1,t2;
	double	sqrt();
	int	i;

	*nout = a = aa = b = bb = ab = 0.0;

	for(i=0;i<nin;i++){
		t1 = arr1[i];
		t2 = arr2[i];
		if(t1 == -1.0 || t2 == -1.0)
			continue;
		(*nout)++;
		a += t1;
		aa += t1*t1;
		b += t2;
		bb += t2*t2;
		ab += t1*t2;
		}

	t1 = ab - a*b/(*nout);
	t2 = (aa - a*a/(*nout))*(bb - b*b/(*nout));
	return(t1/sqrt(t2));
}	
