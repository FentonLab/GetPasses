#include <stdlib.h>
#include <stdio.h>
#include "time_series.h"

void    get_ts_data(TSFILE *tsfile, TIME_SERIES *ts)
{       
	FILE 	*fp;
	unsigned char xb, yb, n_spikes, *byte;
	unsigned int	spike_time;
	int	x, y, i, j, n, n_samps, total_spks = 0, spk_cnt, ReverseBytes;
	char	line[256], string[10], keyword[256], *get_line();
	long	data_offset;
	float	dummy, reduce_x, reduce_y;

	ReverseBytes = tsfile->reverse_bytes;
	reduce_x = tsfile->reduce_x;
	reduce_y = tsfile->reduce_y;

	fp = fopen(tsfile->file, "r");
	if(fp == NULL){
		printf("Can't open %s\n", tsfile->file);
		exit(-1);
	}

	tsfile->samps_per_sec = 0.0;

	if(fscanf(fp,"%d\n", &n) != EOF){
		n--; // to account for the 1st line
		for(i = 0; i < n; i++){	
			get_line(fp, line);
			sscanf(line,"%s%f", &keyword, &dummy);
			if(!strcmp(keyword,"%SAMPLING_INTERVAL(samps/sec)"))
				tsfile->samps_per_sec = (float)dummy;
			if(!strcmp(keyword,"%SCALE_Y(RatioTracktoMapPixels)"))
			 	tsfile->scale_y = dummy;
			if(!strcmp(keyword,"%SCALE_X(RatioTracktoMapPixels)")){
				tsfile->scale_x = dummy;
			}
		}
	}else{
		printf("file is corrupt\n");
		exit(-1);
	}

	if(tsfile->samps_per_sec == 0){
		printf("Did not find %%SAMPLING_INTERVAL(samps/sec) in header. File is corrupt\n");
		exit(-1);
	}
	if(tsfile->scale_y == 0.0){
		printf("Did not find %%SCALE_Y(RatioTracktoMapPixels) in header. File is corrupt\n");
		exit(-1);
	}
	if(tsfile->scale_x == 0.0){
		printf("Did not find %%SCALE_X(RatioTracktoMapPixels) in header. File is corrupt\n");
		exit(-1);
	}

	data_offset = ftell(fp);
	n_samps = 0;
	total_spks = 0;
	while(fread(&xb,sizeof(unsigned char), 1, fp) != 0){
		fread(&yb,sizeof(unsigned char), 1, fp);
		fread(&n_spikes,sizeof(unsigned char), 1, fp);

		//printf("%d\t%d\t%d\t%d\n",x,y,n_spikes, n_samps);

		for(i = 0; i < n_spikes; i++){
			fread(&spike_time,sizeof(unsigned int), 1, fp);
		}
		total_spks += n_spikes;
		n_samps++;
	}
	tsfile->n_samps = n_samps;
	tsfile->Xsz = (int)((float)(tsfile->Xsz) / reduce_x);
	tsfile->Ysz = (int)((float)(tsfile->Ysz) / reduce_y);


	for(y = 0; y <  tsfile->Ysz; y++)
		for(x = 0; x <  tsfile->Xsz; x++)
			ts->fm.time_map[y][x] = ts->fm.spk_map[y][x] = 0;

	ts->n_spks = total_spks;
	ts->p = (POSITION *)calloc((size_t)(n_samps), sizeof(POSITION));
	ts->pos = (POSITION *)calloc((size_t)(n_samps), sizeof(POSITION));
        ts->s = (int *)calloc((size_t)(n_samps), sizeof(int));
        ts->tm = (unsigned long *)calloc((size_t)total_spks, sizeof(unsigned long));
        if((ts->p == NULL) || (ts->pos == NULL) || (ts->s == NULL) || (ts->tm == NULL)){
                fprintf(stderr,"Can't allocate timeseries\n");
                exit(-2);
        }

	fseek(fp, (long) data_offset, SEEK_SET);
	if(ReverseBytes)
		byte = (unsigned char  *)&spike_time;

	spk_cnt = 0;
	for(i=0; i < n_samps; i++){
		fread(&xb,sizeof(unsigned char), 1, fp);
		fread(&yb,sizeof(unsigned char), 1, fp);
		fread(&n_spikes,sizeof(unsigned char), 1, fp);

		ts->pos[i].x = (int)xb;
		ts->pos[i].y = (int)yb;
		ts->s[i] = (int)n_spikes;
		// printf("%d\t%d\t%d\t%d\t",i, ts->p[i].x,ts->pos[i].y,ts->s[i]);

		ts->p[i].x = x = (int)((float)xb / reduce_x);
		ts->p[i].y = y = (int)((float)yb / reduce_y);
                ts->fm.time_map[y][x]++;
                ts->fm.spk_map[y][x] += n_spikes;

		for(j = 0; j < n_spikes; j++){
			if(ReverseBytes){
				fread(byte+3, sizeof(unsigned char), 1, fp);
				fread(byte+2, sizeof(unsigned char), 1, fp);
				fread(byte+1, sizeof(unsigned char), 1, fp);
				fread(byte, sizeof(unsigned char), 1, fp);
			}else{
				fread(&spike_time, sizeof(unsigned int), 1, fp);
			}
			ts->tm[spk_cnt] = spike_time;
			// printf("%d\t", ts->tm[spk_cnt]);
			spk_cnt++;
		}
		// printf("\n");
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
