#include <stdlib.h>
#include <stdio.h>
#include "km.h"
#include "time_series.h"

void    get_km_data(TSFILE *tsfile, TIME_SERIES *ts, float samps_per_sec)
{       
	FILE 	*fp;
	unsigned char x, y, bx, by, n_spikes1;
	int	i, j, n, n_samps, total_spks = 0, spk_cnt; 
	float	reduce_x, reduce_y;
	unsigned char	data[6];
	unsigned int	time_step;


	reduce_x = tsfile->reduce_x;
	reduce_y = tsfile->reduce_y;

	fp = fopen(tsfile->file, "r");
	if(fp == NULL){
		printf("Can't open %s\n", tsfile->file);
		exit(-1);
	}

 	tsfile->scale_y = 1.0;
	tsfile->scale_x = 1.0;
	tsfile->samps_per_sec = samps_per_sec;
	time_step = (int)(10000 / samps_per_sec);

	n_samps = 0;
	total_spks = 0;
	fseek(fp, 0, SEEK_END);
	n_samps = ftell(fp); // total bytes in file
	n_samps /= 6; // 6 bytes per sample

	tsfile->n_samps = n_samps;
	tsfile->Xsz = (int)((float)(tsfile->Xsz) / reduce_x);
	tsfile->Ysz = (int)((float)(tsfile->Ysz) / reduce_y);

	for(y = 0; y <  tsfile->Ysz; y++)
		for(x = 0; x <  tsfile->Xsz; x++)
			ts->fm.time_map[y][x] = ts->fm.spk_map[y][x] = 0;


	ts->p = (POSITION *)calloc((size_t)(n_samps), sizeof(POSITION));
	ts->pos = (POSITION *)calloc((size_t)(n_samps), sizeof(POSITION));
        ts->s = (int *)calloc((size_t)(n_samps), sizeof(int));

        if((ts->p == NULL) || (ts->pos == NULL) || (ts->s == NULL)){
                fprintf(stderr,"Can't allocate timeseries\n");
                exit(-2);
        }

	rewind(fp);

	total_spks = 0;
	for(i=0; i < n_samps; i++){
		fread(data,sizeof(unsigned char), 6, fp);
                y = *(data + BYTE_FOR_RED_Y);
                x = *(data + BYTE_FOR_RED_X);
                by = *(data + BYTE_FOR_BLU_Y);
                bx = *(data + BYTE_FOR_BLU_X);
                n_spikes1 = (int)(*(data + 4) & 15);
               	// n_spikes2 = (int)(*(data + 4) >> 4);
                // n_spikes3 = (int)(*(data + 5) & 15);
               	// n_spikes4 = (int)(*(data + 5) >> 4);

		ts->pos[i].x = (int)x;
		ts->pos[i].y = (int)y;
		ts->s[i] = (int)n_spikes1;
		// printf("%d\t%d\t%d\t%d\t",i, ts->p[i].x,ts->pos[i].y,ts->s[i]);

		ts->p[i].x = x = (int)((float)x / reduce_x);
		ts->p[i].y = y = (int)((float)y / reduce_y);
                ts->fm.time_map[y][x]++;
                ts->fm.spk_map[y][x] += n_spikes1;

		total_spks += n_spikes1;
	}

	// these are set to conform to the ts format
	ts->n_spks = total_spks;
	ts->tm = (unsigned long *)calloc((size_t)total_spks, sizeof(unsigned long));
        if(ts->tm == NULL){
                fprintf(stderr,"Can't allocate timeseries spike times\n");
                exit(-2);
        }

	// spike time info is not in the km format, estimate the spike times for the ts format
	srand48((long)total_spks);
	rewind(fp);
	spk_cnt = 0;
	for(i=0; i < n_samps; i++){
		fread(data,sizeof(unsigned char), 6, fp);
                n_spikes1 = (int)(*(data + 4) & 15);
		for(j=0; j < n_spikes1; j++){
			ts->tm[spk_cnt] = i * time_step + (int)(drand48() * time_step); 
			spk_cnt++;
		}
	}

	for(y = 0; y <  tsfile->Ysz; y++)
		for(x = 0; x <  tsfile->Xsz; x++){
			if(ts->fm.time_map[y][x] > 0){
				ts->fm.rate_map[y][x] = (float)ts->fm.spk_map[y][x] / (float)ts->fm.time_map[y][x] * (float)tsfile->samps_per_sec;
			}else
				ts->fm.rate_map[y][x] =  -1.0;
		}
	return;
}
