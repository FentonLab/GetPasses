#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
// #include <unistd.h>
#include "tsf.h"

void	get_int_map(fp, m, X, Y)
	FILE	*fp;
	INT_MAP	m;
	int	X, Y;
{
	int	x,y;

	for(y = 0; y < Y; y++){
		for(x = 0; x < X; x++){
			if((fscanf(fp, "%d", &(m[y][x])) != 1)){
				fprintf(stderr,"Error reading map element Y= %d X= %d\n",y,x);
				exit(2);
			}
		}
	}
	return;
}
void	get_key(fp, k, n)
	FILE	*fp;
	KEY	*k;
	int	n;
{
	int	i;
	char	s[256], *get_line();

	for(i = n-1; i >= 0; i--){
		if(fscanf(fp,"%f%f", &(k[i].max), &(k[i].medn)) != 2){
			fprintf(stderr,"Error reading key class %d\n",i);
			exit(2);
		}
	}
	(void)get_line(fp, s); 	/* new line */
	return;
}

void	get_float_map(fp, r, X, Y)
	FILE		*fp;
	FLOAT_MAP	r;
	int	X, Y;
{
	int	x,y;

	for(y = 0; y < Y; y++){
		for(x = 0; x < X; x++){
			if((fscanf(fp, "%f", &(r[y][x])) != 1)){
				fprintf(stderr,"Error reading rate element Y= %d X= %d\n",y,x);
				exit(2);
			}
		}
	}
	return;
}

void	minus1_to_0(map, Y, X)
	INT_MAP	map;
	int	Y, X;
{
	int	y, x;

	for(y=0; y < Y;y++)
		for(x=0; x < X;x++)
			if(map[y][x] == -1)
				map[y][x] = 0;
	return;
}

void	make_time_map(map, p, n_samps, Ysz, Xsz)
INT_MAP	map;
POSITION	*p;
int		n_samps;
{
	int	i, x, y;

	for(y = 0; y < Ysz; y++)
		for(x = 0; x < Xsz; x++)
			map[y][x] = 0;

	for(i = 0; i < n_samps; i++){	
		y = p[i].y;
		x = p[i].x;
printf("%d %d\n",y,x);
		if(!y && !x)
			continue;
		map[y][x]++;
	}
	return;
}

void	make_spike_map(map, p, spk, n_samps, Ysz, Xsz)
INT_MAP	map;
POSITION	*p;
int		n_samps, *spk;
{
	int	i, x, y;

	for(y = 0; y < Ysz; y++)
		for(x = 0; x < Xsz; x++)
			map[y][x] = 0;

	for(i = 0; i < n_samps; i++){	
		y = p[i].y;
		x = p[i].x;
		if(!y && !x)
			continue;
		map[y][x] += spk[i];
	}
	return;
}

void	make_rate_map(spk, time, rate, Ysz, Xsz, samps_per_sec)
INT_MAP	spk, time;
FLOAT_MAP	rate;
int	Ysz, Xsz;
float	samps_per_sec;
{
	int	x, y;

	for(y = 0; y < Ysz; y++){
		for(x  = 0; x < Xsz; x++){
			if(time[y][x] > 0){
				rate[y][x] = (float)spk[y][x] / (float)time[y][x] * samps_per_sec;
			}else
				rate[y][x] = -1.0;
		}
	}
	return;
}
void	make_field_map(fields, rate, field, cent_y, cent_x, Ysz, Xsz, Threshold)
	FLOAT_MAP	rate;
	INT_MAP		field, fields;
	float		*cent_y, *cent_x;
	int		Ysz, Xsz;
	float		Threshold;
{
	int		x, y;
	float		wx, wy, w;

	w = wx = wy = 0.0;
	for(y = 0; y < Ysz; y++){
		for(x = 0; x < Xsz; x++){
			if((rate[y][x] > Threshold) && (fields[y][x] == 1)){ /* only consider largest field */
				field[y][x] = 1;
				wy += (float)(rate[y][x] * y);
				wx += (float)(rate[y][x] * x);
				w += rate[y][x];
			}else{
				field[y][x] = 0;
			}
		}
	}
	*cent_y = (wy / w);
	*cent_x = (wx / w);

/*
	for(y = 0; y < Ysz; y++){
		for(x = 0; x < Xsz; x++)
			printf("%d ",field[y][x]);
		printf("\n");
	}
	printf("\n");
*/
	return;
}
void get_avg_rate(avg_rate, spike, time, Ysz, Xsz, samps_per_sec)
	INT_MAP	spike, time;
	int	Ysz, Xsz;
	float	*avg_rate, samps_per_sec;
{
	int	x, y, n, s, t;
	float	rate;

	s = t = n = 0;
	for(y = 0; y < Ysz; y++){
		for(x  = 0; x < Xsz; x++){
			if(time[y][x] > 0){
				s += spike[y][x];		
				t += time[y][x];
				n++;
			}
		}
	}
	*avg_rate = (float)s / (float)t * samps_per_sec;

	return;
}
