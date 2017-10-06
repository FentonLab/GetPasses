#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <math.h>
#include "tsf.h"
#include "dif.h"
#include "get_passes.h"

UNIT *	get_approach_departure_maps(unit, ts, n_samps, h, n_approach_samps, n_departure_samps, t, dwell_criterion, trk_file, target, reduce_x, reduce_y, max_spks, max_psth_bins, samps_per_sec)	
UNIT		*unit;
TIME_SERIES	ts;
HEADER 		h;
int	n_samps, n_approach_samps, n_departure_samps, dwell_criterion, max_spks;
float	samps_per_sec;
TARGET	t;
char	*trk_file, *target;
float	reduce_x, reduce_y;
{
	int	i, j, k, x, y, n_units, n_frames, Ysz, Xsz, direction, n_visits = 0;
	FILE	*afp, *dfp, *zfp;
	static UNIT	*u;
	POSITION	*p;
	void	write_trk(), initialize_maps(), update_maps(), calc_psth(), calc_z(), calc_psdh();
	char	f[256], command[256];
	int	enter_ts, start_ts, stop_ts, ts_offset, total_psth_bins, entrance_samp;
	int	frame, start_samp, stop_samp, approach_start, approach_end, depart_start, depart_end;

	n_units = h.n_units;
	n_frames = h.n_maps;
	Ysz = h.Ysz;
	Xsz = h.Xsz;

	u = (UNIT *)(calloc((size_t)n_units, sizeof(UNIT)));
	if(u == NULL){
		fprintf(stderr,"Can't allocate UNITs in categorize_path\n");
		exit(2);
	}
	total_psth_bins = max_psth_bins * 2;
	ts_offset = max_psth_bins * COUNTS_PER_SAMP;

	for(i=0; i < n_units; i++){
                if((u[i].fm = (FRAME_MAP *)calloc((unsigned)(n_frames * 2), sizeof(FRAME_MAP))) == NULL){
                        fprintf(stderr, "Can't allocate FRAME_MAPs categorize_path\n");
                        exit(2);
                }
		for(j=0; j < (n_frames * 2); j++){
			frame = j / 2;
                        if ((u[i].fm[j].isi = (int *)calloc((unsigned)max_spks, sizeof(int))) == NULL){
                                fprintf(stderr, "Can't allocate %d isi\n", max_spks);
                                exit(2);
                        }
                        u[i].fm[j].isi_n = 0;

			if(!(j % 2)){	/* only one histogram type per ref frame */
                        	if ((u[i].fm[frame].psth = (int *)calloc((unsigned)total_psth_bins, sizeof(int))) == NULL){
                               	 fprintf(stderr, "Can't allocate %d psth bins\n", total_psth_bins);
                               	 exit(2);
                        	}
                        	u[i].fm[frame].psth_n = 0;

                       		 if ((u[i].fm[frame].a_psdh = (int *)calloc((unsigned)(MAX_PSDH_BINS / PSDH_BIN_SZ) + 1, sizeof(int))) == NULL){
                                	fprintf(stderr, "Can't allocate %d psth bins\n", (MAX_PSDH_BINS / PSDH_BIN_SZ));
                          	     	exit(2);
                       		 }
                        	if ((u[i].fm[frame].d_psdh = (int *)calloc((unsigned)(MAX_PSDH_BINS / PSDH_BIN_SZ) + 1, sizeof(int))) == NULL){
                               		fprintf(stderr, "Can't allocate %d psth bins\n", (MAX_PSDH_BINS / PSDH_BIN_SZ));
                                	exit(2);
                        	}
			}
                }
	}
	for(k=0; k < n_units; k++){
		sprintf(f,"%s.U%d.%s.TRK", trk_file, k + 1, target);
		sprintf(command,"cp -f /usr/local/progs/header.TRK %s",f);
		system(command);
		u[k].fm[0].afp = fopen(f, "a");
		if(u[k].fm[0].afp == NULL){
			fprintf(stderr,"Can't open output file: %s\n",f);
			exit(-1);
		}
	}
	write_trk(ts, 0, n_samps, u, n_units, n_frames); 
	for(k=0; k < n_units; k++)
		fclose(u[k].fm[0].afp);

	if(!t.enabled)
		return(u);

	initialize_maps(u, n_units, n_frames, Ysz, Xsz);
	for(i=0; i < n_frames; i++){
		n_visits = 0;
		for(k=0; k < n_units; k++){
			sprintf(f,"%s.U%dA%d.%sa.TRK", trk_file, k + 1, i + 1, target);
			sprintf(command,"cp -f /usr/local/progs/header.TRK %s",f);
			system(command);
			u[k].fm[i].afp = fopen(f, "a");
			if(u[k].fm[i].afp == NULL){
				fprintf(stderr,"Can't open output file: %s\n",f);
				exit(-1);
			}
			sprintf(f,"%s.U%dA%d.%sd.TRK", trk_file, k + 1, i + 1, target);
			sprintf(command,"cp -f /usr/local/progs/header.TRK %s",f);
			system(command);
			u[k].fm[i].dfp = fopen(f, "a");
			if(u[k].fm[i].dfp == NULL){
				fprintf(stderr,"Can't open output file: %s\n",f);
				exit(-1);
			}
		}
		sprintf(f,"%s.A%d.%s.Z", trk_file, i + 1, target);
		zfp = fopen(f, "w");
		if(zfp == NULL){
			fprintf(stderr,"Can't open output file: %s\n",f);
			exit(-1);
		}

		p = ts.rf[i].p;	
		approach_start = approach_end = depart_start = depart_end = 0;
		for(j = n_approach_samps; j < n_samps; j++){
			x = (int)(p[j].x / reduce_x);
			y = (int)(p[j].y / reduce_y);
			if(t.map[y][x]){     /* entered target */	
				k = j;

				while(t.map[y][x] && (k < (j + dwell_criterion)) && (k < n_samps)) {	 
					x = (int)(p[k].x / reduce_x);
					y = (int)(p[k].y / reduce_y);
					k++;
				}

				if(k == (j + dwell_criterion)){	/* dwell criterion met */
					direction = APPROACH; /* entrance */
					
					approach_end = j;
					approach_start = approach_end - n_approach_samps;

					enter_ts = (j * COUNTS_PER_SAMP);
					start_ts = enter_ts - ts_offset; 
					stop_ts = enter_ts + ts_offset; 
					calc_psth(ts.unit, u, n_units, i, start_ts, enter_ts, stop_ts, max_psth_bins);

					/* calc PSDH */
					start_samp = j - 100;
					if(start_samp < 0)
						start_samp = 0;
					stop_samp = j;
					calc_psdh(ts, u, n_units,t, i, start_samp, stop_samp, direction);
					calc_z(p, ts.unit, (j - n_approach_samps), j, unit, n_units, i, direction, zfp, reduce_x, reduce_y, n_visits, samps_per_sec); 

					/* find exit from target */
					j += dwell_criterion;
					entrance_samp = j;

					x = (int)(p[j].x / reduce_x);
					y = (int)(p[j].y / reduce_y);
				 	while((t.map[y][x]) && (j < n_samps)){
						x = (int)(p[j].x / reduce_x);
						y = (int)(p[j].y / reduce_y);
						j++;
					}

					depart_start = j;
					update_maps(p, ts.unit, (entrance_samp - n_approach_samps), entrance_samp, entrance_samp, u, n_units, i, direction, reduce_x, reduce_y); 

					direction = NEUTRAL;
					/* define neutral pass as path from last departure end to current approach start */
					if(depart_end < 0)
						depart_end = 0;
					if(approach_start < 1)
						approach_start = 1;
					calc_z(p, ts.unit, (depart_end + 1), (approach_start - 1), unit, n_units, i, direction, zfp, reduce_x, reduce_y, n_visits, samps_per_sec); 

					/* define neutral pass as path centered between last departure end and this approach start 
					start_samp = depart_end + (int)((approach_start - depart_end) / 2) - (int)(n_approach_samps / 2);
					if(start_samp <= depart_end)
						start_samp = depart_end +1;
					stop_samp = start_samp + (n_approach_samps);
					if(stop_samp >= approach_start)
						stop_samp = approach_start - 1;
					calc_z(p, ts.unit, start_samp, (stop_samp), unit, n_units, i, direction, zfp, reduce_x, reduce_y, n_visits, samps_per_sec); 
					*/

					/* departure */
					direction = DEPARTURE; 
					if((j + n_departure_samps) < n_samps){
						depart_end = j + n_departure_samps;
						update_maps(p, ts.unit, j, (j + n_departure_samps), (j + n_departure_samps), u, n_units, i, direction, reduce_x, reduce_y); 
						calc_z(p, ts.unit, j, (j + n_departure_samps), unit, n_units, i, direction, zfp, reduce_x, reduce_y, n_visits, samps_per_sec); 
					}else{
						depart_end = n_samps - 1;
						update_maps(p, ts.unit, j, n_samps, n_samps, u, n_units, i, direction, reduce_x, reduce_y); 
						calc_z(p, ts.unit, j, n_samps, unit, n_units, i, direction, zfp, reduce_x, reduce_y, n_visits, samps_per_sec); 
					}

					/* calc PSDH */
					start_samp = j;
					stop_samp = j + 100;
					if(stop_samp > n_samps)
						stop_samp = n_samps;
					calc_psdh(ts, u, n_units,t, i, start_samp, stop_samp, direction);

					j += n_departure_samps;
					n_visits++;
				}else{	
					j = k;
				}
			}
		}
		fprintf(stdout,"%s\t Frame: %d\t%d visits\n", trk_file, i+1, n_visits);

		for(k=0; k < n_units; k++){
			fclose(u[k].fm[i].afp);
			fclose(u[k].fm[i].dfp);
		}
		fclose(zfp);
	}
	return(u);
}

void	initialize_maps (u, n_units, n_frames, Ysz, Xsz)
        UNIT            *u;
	int		n_units, n_frames, Ysz, Xsz;
{
	int	i, j, map_num, x, y, n_maps;
	
	/* (Frame 1) appraoch, depart (Frame 2) approach, depart ... */
	n_maps = 2 * n_frames;

	for(i=0; i < n_frames; i++){
		for(j=0; j < n_units; j++){
			for(map_num=0; map_num < n_maps; map_num++){
				for(y = 0; y < Ysz; y++){
					for(x = 0; x < Xsz; x++){
						u[j].fm[map_num].time_map[y][x] = 0;		
						u[j].fm[map_num].spk_map[y][x] = 0;
					}
				}
			}
		}
	} 
	return;
}

void	calc_z(p, unit, start, stop, u, n_units, frame, direction, fp, reduce_x, reduce_y, pass_number, samps_per_sec) 
        POSITION   *p;
        UNIT_LIST  unit[]; /* this contains the time series info for each unit */
	int	start, stop, n_units, frame, direction, pass_number;
	float	 samps_per_sec;
        UNIT	*u;	/* this is the structure containing the ref frame info for the whole session */
	FILE	*fp;
	float	reduce_x, reduce_y;
{
	int	i, j, k, x0, y0, x, y; 
	int	samps, distance_dt = DEFAULT_DISTANCE_dT - 1;
	float	obs_spks, exp_spks, z;
	double	sqrt(), hypot(), distance;


	for(k=0; k < n_units; k++){
		obs_spks = exp_spks = 0.0;
		distance = 0.0;
		samps = stop - start;
		if(samps < 0)
			samps = 0;

		switch(direction){
		case APPROACH:
			fprintf(fp, "%s\ta%d\t%d\t", u[k].fm[frame].file, pass_number, samps);
			break;
		case DEPARTURE:
			fprintf(fp, "%s\td%d\t%d\t", u[k].fm[frame].file, pass_number, samps);
			break;
		case NEUTRAL:
			fprintf(fp, "%s\tn%d\t%d\t", u[k].fm[frame].file, pass_number, samps);
			break;
		default:
			fprintf(stderr,"Error in calc_z. Undefined 'direction'\n");
			exit(0);
		}

		x0 = x = (int)(p[start].x / reduce_x);
		y0 = y = (int)(p[start].y / reduce_y);
		for(j=start, samps = 0; j < stop; j++, samps++){

			x = (int)(p[j].x / reduce_x);
			y = (int)(p[j].y / reduce_y);

			/* shouldn't be necessary but it is because Kaminsky either set what was "not in the arena to 0 in the time map or the conversion of the track coordinates from the DIF file is not accurate ... sigh */
			if(!u[k].fm[frame].time_map[y][x]){
				continue;
			}

			/* shouldn't be necessary but what the hell */
			if(!x && !y){
				continue;
			}
		
			obs_spks += (float)unit[k].s[j];		
			exp_spks +=(float) (u[k].fm[frame].rate_map[y][x]) / samps_per_sec;		

			if((samps % DEFAULT_DISTANCE_dT) == distance_dt){
				distance += hypot((double)(x - x0), (double)(y - y0));
				x0 = x;
				y0 = y;
			}
		}
		if(exp_spks > 0.0)
			z = (obs_spks - exp_spks) / (float)sqrt(exp_spks);
		else
			z = 0.0;
		fprintf(fp, "%3.2lf\t%3.0f\t%3.2f\t%3.2f\n", distance, obs_spks, exp_spks, z);
	}
	return;
}

void	write_trk(ts, start, stop, u, n_units, n_frames)
        TIME_SERIES ts;
	int	start, stop, n_units, n_frames;
        UNIT	*u;
{
	int	i, j, k;

	for(j=start; j < stop; j++){
		for(k=0; k < n_units; k++){
			for(i = 0; i < n_frames; i ++){
				fprintf(u[k].fm[0].afp, "%d\t%d\t",ts.rf[i].p[j].x, ts.rf[i].p[j].y);
				if(n_frames == 1)	/* just to make a standardized 2 frame trk file */
					fprintf(u[k].fm[0].afp, "%d\t%d\t",0, 0);
			}
			fprintf(u[k].fm[0].afp, "%d\n",ts.unit[k].s[j]);
		}
	}
	return;
}

void	update_maps(p, unit, start, enter, stop, u, n_units, frame, direction, reduce_x, reduce_y)
        POSITION   *p;
        UNIT_LIST  unit[];
	int	start, enter, stop, n_units, frame, direction;
        UNIT	*u;
	float	reduce_x, reduce_y;
{
	int	i, j, k, x, y, map_num;
	void	get_isi();

	for(j=start; j < stop; j++){
		/* don't count samples in the target area */
		if(j > enter)
			continue;

		x = (int)(p[j].x / reduce_x);
		y = (int)(p[j].y / reduce_y);

		for(k=0; k < n_units; k++){
			if(!frame){ /* room */
				if(!direction){ /* approach */
					fprintf(u[k].fm[frame].afp, "%d\t%d\t%d\t%d\t",p[j].x, p[j].y, 0, 0);
				}else{ 	/*departure */
					fprintf(u[k].fm[frame].dfp, "%d\t%d\t%d\t%d\t",p[j].x, p[j].y, 0, 0);
				}
			}else{	/* arena */
				if(!direction){ /* approach */
					fprintf(u[k].fm[frame].afp, "%d\t%d\t%d\t%d\t",0, 0, p[j].x, p[j].y);
				}else{		/* departure */
					fprintf(u[k].fm[frame].dfp, "%d\t%d\t%d\t%d\t",0, 0, p[j].x, p[j].y);
				
				}
			}

			/*
			if(j == enter){
				if(!direction){ * approach *
					fprintf(u[k].fm[frame].afp, "0\n");
				}else{		* departure *
					fprintf(u[k].fm[frame].dfp, "0\n");
				}
				continue;
			}
			if(j > enter)
				continue;
			*/

			/* (Frame 1) approach, depart (Frame 2) approach, depart ... */
			map_num = frame * 2 + direction;
			(u[k].fm[map_num].time_map[y][x])++;		
			(u[k].fm[map_num].spk_map[y][x]) += unit[k].s[j];		

			if(!direction)	/* approach */
				fprintf(u[k].fm[frame].afp, "%d\n", unit[k].s[j]);
			else 		/* departure */
				fprintf(u[k].fm[frame].dfp, "%d\n", unit[k].s[j]);

			if(unit[k].s[j])
				get_isi(u[k].fm[map_num].isi, &(u[k].fm[map_num].isi_n), unit[k].ts, unit[k].n, j);
		}
	}
	/* make a marker between the passes */
	for(k=0; k < n_units; k++){
		if(!direction){ /* approach */
			fprintf(u[k].fm[frame].afp, "%d\t%d\t%d\t%d\t%d\n",0, 0, 0, 0, 0);
		}else{
			fprintf(u[k].fm[frame].dfp, "%d\t%d\t%d\t%d\t%d\n",0, 0, 0, 0, 0);
		}
	}
	return;
}

void	make_target_map(h, t)
	HEADER h;
	TARGET	*t;
{
	double hypot();
	int x, y, Xsz, Ysz;

	Xsz = h.Xsz;
	Ysz = h.Ysz;
	for(y = 0; y < Ysz; y++){
		for(x = 0; x < Xsz; x++){
			if(hypot(x - t->x, y - t->y) > t->r)
				t->map[y][x] = 0;
			else{
				t->map[y][x] = 1;
			}
		}
	}
	return;
}

void	write_maps(h, u, unit, n_units, n_frames, target, s_file, t_file, r_file, i_file, psth_file, psdh_file, total_psth_bins, samps_per_sec, target_enabled)
	HEADER	h;
	int	n_units, n_frames, total_psth_bins;
        UNIT	*u, *unit;
	char	*target, *s_file, *t_file, *r_file, *i_file, *psth_file, *psdh_file;
	float	samps_per_sec;
{
	FILE	*sfp, *tfp, *rfp, *ifp, *psthfp, *psdhfp;
	char	f[256];
	int	i, j, k, x, y, frame, s, t, n_maps, bin, max_bin, psdh_bins, n_spks;
	float	r;

	if(!target_enabled){
		for(i=0; i < n_units; i++){
			for(frame=0; frame < n_frames; frame++){
				sprintf(f, "%s.U%dA%d", s_file, i+1, frame + 1);
				sfp = fopen(f, "w");
				if(sfp == NULL){
					fprintf(stderr, "Can't open file %s\n", f);
					exit(-1);
				}

				sprintf(f, "%s.U%dA%d", t_file, i+1, frame + 1);
				tfp = fopen(f, "w");
				if(tfp == NULL){
					fprintf(stderr, "Can't open file %s\n", f);
					exit(-1);
				}

				sprintf(f, "%s.U%dA%d", r_file, i+1, frame + 1);
				rfp = fopen(f, "w");
				if(rfp == NULL){
					fprintf(stderr, "Can't open file %s\n", f);
					exit(-1);
				}
				for(y = 0; y < h.Ysz; y++){
					for(x = 0; x < h.Xsz; x++){
						t = unit[i].fm[frame].time_map[y][x];		
						if(t){
							s = unit[i].fm[frame].spk_map[y][x];		
							r = (float)(s) / (float)t * samps_per_sec;
						}else{
							t = s = -1;
							r = -1.0;
						}
						fprintf(sfp, "%d\n", s);
						fprintf(tfp, "%d\n", t);
						fprintf(rfp, "%0.3f\n", r);
					}
				}
			}
		}
		return;
	}


	n_maps = n_frames * 2;
	for(i=0; i < n_units; i++){
		for(j = 0; j < n_maps; j++){
			/* (Frame 1) approach, depart (Frame 2) approach, depart ... */
			frame = (int)(j / 2);	/* two maps (approach and departure) per frame */

			sprintf(f, "%s.U%dA%d.%s", psth_file, i+1, frame + 1, target);
			psthfp = fopen(f, "w");
			if(psthfp == NULL){
				fprintf(stderr, "Can't open file %s\n", f);
				exit(-1);
			}

			if(!(j % 2)){ /* an approach map */ 
				sprintf(f, "%s.U%dA%da.%s", s_file, i+1, frame + 1, target);
				sfp = fopen(f, "w");
				if(sfp == NULL){
					fprintf(stderr, "Can't open file %s\n", f);
					exit(-1);
				}

				sprintf(f, "%s.U%dA%da.%s", t_file, i+1, frame + 1, target);
				tfp = fopen(f, "w");
				if(tfp == NULL){
					fprintf(stderr, "Can't open file %s\n", f);
					exit(-1);
				}

				sprintf(f, "%s.U%dA%da.%s", r_file, i+1, frame + 1, target);
				rfp = fopen(f, "w");
				if(rfp == NULL){
					fprintf(stderr, "Can't open file %s\n", f);
					exit(-1);
				}

				sprintf(f, "%s.U%dA%da.%s", i_file, i+1, frame + 1, target);
				ifp = fopen(f, "w");
				if(ifp == NULL){
					fprintf(stderr, "Can't open file %s\n", f);
					exit(-1);
				}

				sprintf(f, "%s.U%dA%da.%s", psdh_file, i+1, frame + 1, target);
				psdhfp = fopen(f, "w");
				if(psdhfp == NULL){
					fprintf(stderr, "Can't open file %s\n", f);
					exit(-1);
				}
			}else{	 /* a departure map */
				sprintf(f, "%s.U%dA%dd.%s", s_file, i+1, frame + 1, target);
				sfp = fopen(f, "w");
				if(sfp == NULL){
					fprintf(stderr, "Can't open file %s\n", f);
					exit(-1);
				}

				sprintf(f, "%s.U%dA%dd.%s", t_file, i+1, frame + 1, target);
				tfp = fopen(f, "w");
				if(tfp == NULL){
					fprintf(stderr, "Can't open file %s\n", f);
					exit(-1);
				}

				sprintf(f, "%s.U%dA%dd.%s", r_file, i+1, frame + 1, target);
				rfp = fopen(f, "w");
				if(rfp == NULL){
					fprintf(stderr, "Can't open file %s\n", f);
					exit(-1);
				}

				sprintf(f, "%s.U%dA%dd.%s", i_file, i+1, frame + 1, target);
				ifp = fopen(f, "w");
				if(ifp == NULL){
					fprintf(stderr, "Can't open file %s\n", f);
					exit(-1);
				}

				sprintf(f, "%s.U%dA%dd.%s", psdh_file, i+1, frame + 1, target);
				psdhfp = fopen(f, "w");
				if(psdhfp == NULL){
					fprintf(stderr, "Can't open file %s\n", f);
					exit(-1);
				}
			}
					
			for(y = 0; y < h.Ysz; y++){
				for(x = 0; x < h.Xsz; x++){
					t = u[i].fm[j].time_map[y][x];		
					if(t){
						s = u[i].fm[j].spk_map[y][x];		
						r = (float)(s / t * samps_per_sec);
					}else{
						t = s = -1;
						r = -1.0;
					}
					fprintf(sfp, "%d\n", s);
					fprintf(tfp, "%d\n", t);
					fprintf(rfp, "%0.3f\n", r);
				}
			}
			/* write the ISI */
			for(k = 0; k < u[i].fm[j].isi_n; k++)
				fprintf(ifp, "%ld\n",u[i].fm[j].isi[k]);

			if(!(j % 2)){	/* only 1 PSTH per frame */
				/* write the psth */
				bin =  0 - (total_psth_bins / 2);
				max_bin = total_psth_bins;
				for(k = 0;k <= max_bin; k++, bin++){
					fprintf(psthfp, "%d\t%d\n", bin, (int)u[i].fm[frame].psth[k]);
				}
			}

			/* write the psdh */
			psdh_bins = MAX_PSDH_BINS / PSDH_BIN_SZ + 1;
			if(j % 2){ /* departures */
				for(bin = 0, n_spks = 0;bin < psdh_bins; bin++)
					n_spks += u[i].fm[frame].d_psdh[bin]++;
				for(bin = 0;bin < psdh_bins; bin++){
					fprintf(psdhfp, "%d\t%d\t%0.4f\n", bin, u[i].fm[frame].d_psdh[bin], (float)u[i].fm[frame].d_psdh[bin] / (float)n_spks);
				}
			}else{	/* approaches */
				for(bin = 0, n_spks = 0;bin < psdh_bins; bin++)
					n_spks += u[i].fm[frame].a_psdh[bin]++;
				for(bin = psdh_bins;bin--;){
					fprintf(psdhfp, "%d\t%d\t%0.4f\n", 0 - bin, u[i].fm[frame].a_psdh[bin], (float)u[i].fm[frame].a_psdh[bin] / (float)n_spks);
				}
			}

			fclose(sfp);
			fclose(tfp);
			fclose(rfp);
			fclose(ifp);
			fclose(psthfp);
			fclose(psdhfp);
		}
	}

	return;
}

void	get_isi(isi_a, n, ts, n_spks, samp)
	int	*isi_a, *n, n_spks, samp;
	unsigned long	*ts;
{
	int	i;
	int	isi;
	i= 0; 

	while ((i < n_spks) && ((int)(ts[i] / COUNTS_PER_SAMP) < samp)){
		i++; 
	}
			
	while ((i < n_spks) && ((int)(ts[i] / COUNTS_PER_SAMP) < (samp + 1))){
		isi = (ts[i] - ts[i -1]) / COUNTS_PER_MS;	/* isi in ms */
		if(isi < MAX_ISI){
			isi_a[*n] = isi;
			(*n)++;	
		}
		i++; 
	}
	return;
}

void	calc_psth(unit, u, n_units, frame, start, enter, stop, max_bin)
	UNIT_LIST	*unit;
	UNIT		*u;	
	int		n_units, max_bin;
	int		frame, start, enter, stop;
{
	int	i, j, bin;

	for( i = 0; i < n_units; i++){
		/* find first spike */
		j = 0;

		while((unit[i].ts[j] < start) && (j < unit[i].n))
			j++;		
		/* fill in psth */
		while((unit[i].ts[j] < stop) && (j < unit[i].n)){
 			bin = max_bin + (((int)unit[i].ts[j] - enter) / COUNTS_PER_SAMP);

			u[i].fm[frame].psth[bin]++;
			u[i].fm[frame].psth_n++;
			j++;
		}
	}
	return;			
}

void	calc_psdh(ts, u, n_units,t, frame, start, stop, direction)
	TIME_SERIES	ts;
	UNIT		*u;	
	TARGET		t;
	int		n_units;
	int		frame, start, stop, direction;
{
	int	i, j, bin, d;

	for( i = 0; i < n_units; i++){
		for(j = start; j < stop; j++){
			if(ts.unit[i].s[j]){	/* if there is a spike */
 				bin = (int)(hypot(ts.rf[frame].p[j].x - t.x, ts.rf[frame].p[j].y - t.y) + 0.5);
				bin /= PSDH_BIN_SZ;

				if(direction)	/* departure */
					u[i].fm[frame].d_psdh[bin]++;
				else
					u[i].fm[frame].a_psdh[bin]++;
			}
		}
	}
	return;			
}
