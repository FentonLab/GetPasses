#include	"tsf.h"

int	get_cor_arr(cor_arr,t,s,Ysz,Xsz,samp_per_sec)
float	cor_arr[][2], samp_per_sec;
int	Ysz,Xsz;
INT_MAP	s,t;

{
	float	near_val();
	int	i,j,n=0;
	int	spks,time;

	for(i=1;i<Ysz - 1; i++)
		for(j=1;j<Xsz - 1; j++){
			if(t[i][j] == 0);
			else{
				spks = s[i][j];
				time = t[i][j];
				cor_arr[n][0] =(((float) spks)/((float) time))* samp_per_sec;
				cor_arr[n++][1] = near_val(t,s,i,j,samp_per_sec);
				if(cor_arr[n-1][1] == -1.0)
					n--;
				}
			}
	return(n);
}

float	near_val(t,s,y,x,samp_per_sec)
INT_MAP	t;
INT_MAP	s;
int	y,x;
float	samp_per_sec;

{
	int	i,j;
	int	ssum,tsum;

	ssum = tsum = 0;

	for(i=y-1;i<y+2;i++)
		for(j=x-1;j<x+2;j++){
			if(i==y && j==x);
			else{
				ssum += s[i][j];
				tsum += t[i][j];
				}
			}
	if(tsum != 0.0)
		return(((double) ssum)/((double) tsum)*samp_per_sec);
	else
		return(-1.0);
}
