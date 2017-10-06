#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	"time_series.h"
#include	"get_passes.h"

extern	int	optind;
extern	char	*optarg;

main(int argc, char *argv[])
{
	char	*strcpy(), *getenv(), *trunc_str(); 
	void	*calloc(); 
	FILE	*fopen();
	float	get_stat_threshold();
	float	grand_rate;
	int	getopt(); 
	void	instruct(); 
	void	get_ts_data(), get_km_data(), make_spike_map(), make_rate_map(), make_time_map();
	void	find_field_passes(), define_episode_passes();

	TSFILE	tsfile;
	TIME_SERIES	ts;

	INT_MAP	field_map;
	POSITION	*pos;
	int		*spks;	

	TARGET	target;
	FILE	*ifp;
	float	samps_per_sec, reduce_x = REDUCE_X, reduce_y = REDUCE_Y;
	float	ThresholdForCentroid;
	double	targ_ang;
	int	Ysz, Xsz, SamplesPerEpisode = 0;
	int	ReverseBytes = 0, min_samples = -1, n_episodes = 0, n_samps, min_pass_duration;
	int	c,i,j,x,y, dummy; 
	int	max_spks = MAX_SPIKES;
	int	AnalyseEpisodes = 0, calculate_maps = 0;
	int	spot_sz = 2, n_fields = 0;
	int	ts_data_format = 1, print_rate_map = 0;
	str256	file, base_name, episodes_file;
	str256	dat_dir, targ, path_dir;
	EPISODE	*episodes;

	*episodes_file = '\0';
	
	while((c = getopt(argc,argv,"f:k:m:ort:x:y:")) != EOF)
		switch(c){
			case 'f': AnalyseEpisodes = 1;
				strcpy(episodes_file, optarg);
				break;

			case 'k': ts_data_format = 0;
				samps_per_sec = atof(optarg);
				break;

			case 'm': min_samples = atoi(optarg);
				break;

			case 'o': print_rate_map = 1;
				break;

			case 'r': ReverseBytes = 1;
				break;

			case 't': AnalyseEpisodes = 1;
				SamplesPerEpisode = atoi(optarg);
				break;

			case 'x': reduce_x = atof(optarg);
				break;

			case 'y': reduce_y = atof(optarg);
				break;

			default: 	instruct(argv);
		}
	if(argc - optind != 1)
		instruct(argv);
	
	if(getenv("TS_DIR") == NULL){
		(void)fprintf(stderr,"TS_DIR not defined\n");
		SHOW_DIRS
		exit(1);
	}
	strcpy(dat_dir,getenv("TS_DIR"));

	if(getenv("PATH_DIR") == NULL){
		(void)fprintf(stderr,"PATH_DIR not defined. Needed for path data\n");
		SHOW_DIRS
		instruct(argv);
	}
	strcpy(path_dir,getenv("PATH_DIR"));


	strcpy(file,argv[optind]);
	strcpy(base_name, file);
	// trunc_str(base_name,'.');
	strcpy(tsfile.file_name, base_name);
	strcpy(tsfile.path, dat_dir);
	sprintf(tsfile.file,"%s/%s",dat_dir,file);
	tsfile.Xsz = DEFAULT_XSZ;
	tsfile.Ysz = DEFAULT_YSZ;
	tsfile.reduce_x = reduce_x;
	tsfile.reduce_y = reduce_y;
	tsfile.reverse_bytes = ReverseBytes;

	if(ts_data_format)
		get_ts_data(&tsfile, &ts); // makes maps in reduced scale
	else
		get_km_data(&tsfile, &ts, samps_per_sec);// makes maps in reduced scale

	samps_per_sec = tsfile.samps_per_sec;
	if(min_samples == -1)
		min_samples = (int)(samps_per_sec / 5); // 200 ms sampling minimum

	n_samps = tsfile.n_samps;
	Xsz = tsfile.Xsz;
	Ysz = tsfile.Ysz;
	pos = ts.p;
	spks = ts.s;
	grand_rate = (float)ts.n_spks / (float)tsfile.n_samps * tsfile.samps_per_sec;
	
	if(print_rate_map){
		for(y = 0; y <  Ysz; y++){
                	for(x = 0; x <  Xsz; x++){
                               	 fprintf(stdout,"%0.2f\n",ts.fm.rate_map[y][x]);
			}
                }
	}

	// printf("%d %d %f\n", ts.fm.time_map[10][10], ts.fm.spk_map[10][10], ts.fm.rate_map[10][10]);

	if(*episodes_file != '\0'){
		min_pass_duration = 0;
		ifp = fopen(episodes_file, "r");
		if(ifp == NULL){
			fprintf(stderr, "Can't open episodes file\n", episodes_file);
			exit(-1);
		}
		n_episodes = 0;
		while(fscanf(ifp,"%d%d", &dummy, &dummy) != EOF)
			n_episodes++;
		rewind(ifp);

		episodes = (EPISODE *)calloc((size_t)n_episodes, sizeof(EPISODE));
                if(episodes == NULL){
                        fprintf(stderr,"Can't allocate %d EPISODES\n", n_episodes);
                        exit(-1);
                }
		for(i = 0; i < n_episodes; i++){
			fscanf(ifp,"%d%d", &(episodes[i].start), &(episodes[i].stop));
		}
		fclose(ifp);
	}else{
		min_pass_duration = 0;
		n_episodes = (int)(n_samps / SamplesPerEpisode);
		episodes = (EPISODE *)calloc((size_t)n_episodes, sizeof(EPISODE));
		if(episodes == NULL){
			fprintf(stderr,"Can't allocate %d EPISODES\n", n_episodes);
			exit(-1);
		}
		for(i=0; i < n_episodes; i++){
			episodes[i].start = SamplesPerEpisode * i;
			episodes[i].stop = episodes[i].start + SamplesPerEpisode;
		}
	}
	define_episode_passes(spks, pos, n_samps, ts.fm.rate_map, ts.fm.time_map, samps_per_sec, min_samples, min_pass_duration, Ysz, Xsz, path_dir, tsfile.file_name, n_episodes, episodes, grand_rate);  

	return(0);
}

void instruct(argv)
	char	**argv;
{
	fprintf(stderr,"%c\n%s: call with options and a 'ts' (timeseries) file.\n",7,argv[0]);
	fprintf(stderr,"\tInput TS_DIR: .ts data_files\n");
	fprintf(stderr,"\n\tOutput:\n");
	fprintf(stderr,"PATH_DIR: .T\tformat: header: yin xin ycent xcent yout xout\n");
	fprintf(stderr,"           \t for each samp: y x spk rate_in_samp\n");
	fprintf(stderr,"PATH_DIR: .Tz\tformat: obs exp z avg_rate start_sample duration path_length(pix) speed linearity\n");

	fprintf(stderr,"\nOptions:\n");
	fprintf(stderr,"-f\tDefine passes by episodes (start/stop samples) given in the specified file in the pwd\n");
	fprintf(stderr,"-k\tData are in KM format. Provide the value for samps_per_sec.\n");
	fprintf(stderr,"-m\tRequire that a position is sampled a minimum number of samples for estimating rate)\n");
	fprintf(stderr,"  \tA 200 ms minimum is imposed as a default)\n");
	fprintf(stderr,"-o\tPrint rate map values to stdout.\n");
	fprintf(stderr,"-r\tReverse the bytes for interpreting the time stamps.\n\tDo this if the TS file was created on a BIG/LITTLE_ENDIAN machine and is being processed on a machine with the opposite architecture.\n");
	fprintf(stderr,"-t\tDefine passes by episodes (start/stop samples) of a constant length (in sample units)\n");
	fprintf(stderr,"-x\tReduce scale of the raw x coordinates for scaling maps. Default = %0.2f\n", REDUCE_X);
	fprintf(stderr,"-y\tReduce scale of the raw y coordinates  for scaling maps. Default = %0.2f\n", REDUCE_X);
	exit(0);
}

void    find_field_passes(spk, pos, n_samps, rate, field, time, samps_per_sec, centroid_y, centroid_x, Ysz, Xsz, pathdir, base_name, n_episodes,min_samples, spot_sz, grand_rate)
	float		samps_per_sec, centroid_y, centroid_x, grand_rate;
	int		*spk;
	POSITION	*pos;
	FLOAT_MAP	rate;
	INT_MAP		field, time;
        char	*pathdir, *base_name;
        int     n_samps, Ysz, Xsz, n_episodes,min_samples, spot_sz;
{
        FILE    *tfp,*zfp;
        int     i,j, cent_y, cent_x;
        char	out_file[256],temp_name[256];
        INT_MAP hotspot;
        int     x_start,y_start,x_stop,y_stop;

        char    *strcat(),*strcpy();
        FILE    *fopen();
        void    find_paths();

	cent_y = (int)(centroid_y + 0.5);
	cent_x = (int)(centroid_x + 0.5);
        (void)sprintf(out_file,"%s/%s.T",pathdir, base_name);
        if((tfp = fopen(out_file,"w")) == NULL){
                (void)fprintf(stderr,"Can't open output file %s\n",out_file);
                exit(1);
        }
	(void)sprintf(temp_name,"%s/%s.Tz",pathdir,base_name);
        if((zfp = fopen(temp_name,"w")) == NULL){
                (void)fprintf(stderr,"Can't open output file %s\n",temp_name);
                exit(1);
        }

        y_start = cent_y - spot_sz + 1;
        y_stop = cent_y + spot_sz - 1;
        x_start = cent_x - spot_sz + 1;
        x_stop = cent_x + spot_sz - 1;

        for(i= 0;i < Ysz;i++){
                for(j= 0;j < Xsz;j++){
			if((i >= y_start) && (i <= y_stop) && (j >= x_start) && (j <= x_stop)){
                               	hotspot[i][j] = 1;
			}else{
                               	hotspot[i][j] = 0;
			}
		}
	}
	find_paths(tfp, zfp, spk, pos, n_samps, field, time, rate, samps_per_sec, hotspot, min_samples, grand_rate);

        fclose(tfp);
        fclose(zfp);
        return;
}
int     enter_field(pos, field, total_samps, i,y0,x0,y1,x1)
	POSITION	*pos;
	INT_MAP		field;
        int     *i, total_samps;
        unsigned int    *y0, *x0, *y1, *x1;
{
        unsigned int     temp,y=0,x=0,b;
        int     s, n_samps = 0, samps_before = 0, samps_after = 0;

        s = temp = *i;
	x = y = n_samps = 0;

        do{
                *i = s;
                *y0 = (int)((float)pos[s].y);
                *x0 = (int)((float)pos[s].x);
                n_samps++;
                s--;
		if(s < 0)
			break;
                y = (int)((float)pos[s].y);
                x = (int)((float)pos[s].x);
        }while((field[y][x] == 1) || !(y && x));
        s = temp;
        do{
                *y1 = (int)((float)pos[s].y);
                *x1 = (int)((float)pos[s].x);
                s++;
		if(s >= total_samps)
			break;
                n_samps++;
                y = (int)((float)pos[s].y);
                x = (int)((float)pos[s].x);
        }while((field[y][x] == 1) || !(y && x));
        if((*i -= samps_before) < 0)
                *i = 0;
        return(n_samps + samps_before + samps_after -1);
}

void    find_paths(tfp, zfp, spk, pos, n_samps, field, time, rate, samps_per_sec, hotspot, min_samples, grand_rate)
	int	*spk;
	POSITION	*pos;
	INT_MAP		field, time, hotspot;
	FLOAT_MAP	rate;
        int     n_samps, min_samples;
        FILE    *tfp,*zfp;
	float	samps_per_sec, grand_rate;

{
        int     i,b,s;
        unsigned int    x,y;
        unsigned int    x0,y0,x1,y1;
        unsigned int    tx,ty;
        int     in_samp,out_samp, good_samps, d_samps, bin, pass_duration;
	int	min_pass_duration, obs;
	double	z;
	void	calc_z(), calc_locomotion();
	float	speed, linearity, path_length, exp;

	min_pass_duration = SAMPS_PER_SEC;

        for(i=0;i < n_samps;i++){
		y = (int)(pos[i].y);
		x = (int)(pos[i].x);
                if(hotspot[y][x]){
                        pass_duration = enter_field(pos, field, n_samps, &i,&y0,&x0,&y1,&x1);
                        in_samp = i;
                        out_samp = in_samp + pass_duration - 1;
                        if((i + (pass_duration)) > n_samps)
                                break;
                        if(pass_duration < min_pass_duration)
                                break;
                        if(check_sampling(pos, in_samp, out_samp, &good_samps, time, min_samples))
                                break;

                        fprintf(tfp,"%d\t%d  %d\t%d  %d\t%d  %d\n",good_samps,y0,x0,y,x,y1,x1);

			calc_z(i, pass_duration, pos, spk, rate, &obs, &exp, &z, samps_per_sec);
			calc_locomotion(i, pass_duration, pos, samps_per_sec, &path_length, &speed, &linearity);

                        (void)fprintf(zfp,"%d\t%10.8f\t%10.8lf\t%8.5f\t%d\t%8.3f\t%8.3f\t%8.3f\n",obs, exp, z, grand_rate, pass_duration, path_length, speed, linearity);

                        for(s=0;s < pass_duration;s++,i++){
				ty = (int)(pos[i].y);
				tx = (int)(pos[i].x);
                                if(ty && tx){
                                        (void)fprintf(tfp,"%d\t%d\t%d\t%2.6f\n",ty,tx,spk[i],rate[ty][tx]/ (double)samps_per_sec);
                                }
                        }

                }
        }
        return;
}

int     check_sampling(pos, in_samp, out_samp, good_samps,time, min_samples)
	POSITION *pos;
        int     in_samp, out_samp, *good_samps, min_samples;
        INT_MAP	time;
{
        int     i, b, y, x, start, stop;

        *good_samps = 0;
        for(i= in_samp;i <= out_samp;i++){
                y = (int)(pos[i].y);
                x = (int)(pos[i].x);
                if(y && x){
                        (*good_samps)++;
                        if(time[y][x] < min_samples)
                                return(1);
                }
        }
        return 0;
}

/*
int     test_run_speed(start, stop, integ_time, min_distance)
        int     start, stop, integ_time;
        double  min_distance;
{
        int     i, b, y, x, y0 = 0, x0 = 0;
        double hypot(), dist;
        unsigned int    byte[BYTES_PER_SAMP_2B];

        min_distance /= 4.0;
        for(i=start;i <= stop;i += (BYTES_PER_SAMP_2B * integ_time)){
                for(b=0;b<BYTES_PER_SAMP_2B;b++)
                        byte[b] = sess.data[i + b];

                if(reduce_flag){
                        y = ((byte[0] & 63) >> 1);
                        x = ((byte[1] & 63) >> 1);
                        if(y && x){
                                y += 16;
                                x += 16;
                        }
                }else{
                        y = (byte[0] & 63);
                        x = (byte[1] & 63);
                }

                if(i > start){
                        dist = hypot((double)(x0 - x),(double)(y0 - y));
                        if(dist < min_distance)
                                return(0);
                }
                x0 = x;
                y0 = y;
        }
        return(1);
}
*/

void	calc_z(start, n_samps, pos, spk, rate, obs, exp, z, samps_per_sec)
	int	start, n_samps;
	POSITION	*pos;
	int		*spk, *obs;
	FLOAT_MAP	rate;
	double 		*z;
	float		*exp, samps_per_sec;
{
	int 	x, y, i, stop;
	
	*obs = *exp = 0.0;
	stop = start + n_samps;
	for(i = start; i < stop; i++){
		y = (int)(pos[i].y);
		x = (int)(pos[i].x);
		*obs += spk[i];
		*exp += (float)(rate[y][x] / samps_per_sec);	
// printf("XX%f\t%f\t%f\n", *exp, rate[y][x], samps_per_sec);
	}
	*z = ((double)*obs - (double)*exp) / sqrt((double)*exp);

	return;
}
void	calc_locomotion(start, n_samps, pos, samps_per_sec, path_length, speed, linearity)
	int	start, n_samps;
	POSITION	*pos;
	float		samps_per_sec;
	float		*path_length, *speed, *linearity;
{
	int 	n, i, j, k, stop, dist_interval, x0, y0, y, x;
	double	hypot(), Idist, Ldist;
	
	n = 0;
	dist_interval = (int)(samps_per_sec / 2.0);
	*path_length = *linearity = *speed = 0.0;
	Idist = 0.0;
	stop = start + n_samps;

	for(j = start, i = start + 1; i < stop; i++, j++){
		y = (int)(pos[i].y);
		x = (int)(pos[i].x);
		if(y && x){
			k=j;
			while(k >=start){
				y0 = (int)(pos[k].y);
				x0 = (int)(pos[k].x);
				if(x0 && y0){
					Idist += hypot((double)(y - y0), (double)(x - x0)); 		
					break;
				}
				k--;
			}
			if(!((i - start) % dist_interval)){
				k=i;
				while(k >=start){
					if((k - dist_interval) < 0)
						break;

					y0 = (int)(pos[k - dist_interval].y);
					x0 = (int)(pos[k - dist_interval].x);
					if(x0 && y0){
						Ldist = hypot((double)(y - y0), (double)(x - x0)); 		
						if(Idist > 0.0)
							*linearity += (float)(Ldist / Idist);
						*path_length += (float)Ldist;
						n++;
						break;
					}
					k--;
				}
				Idist = 0.0;
			}
		}
	}
		
	*linearity /= n;
	*speed = *path_length / (n * dist_interval) * samps_per_sec;	
	
	return;
}

void    define_episode_passes(spk, pos, n_samps, rate, time, samps_per_sec, min_samples, min_pass_duration, Ysz, Xsz, pathdir, base_name, n_episodes, episodes, grand_rate)
	float		samps_per_sec, grand_rate;
	int		*spk;
	POSITION	*pos;
	FLOAT_MAP	rate;
	INT_MAP		time;
        char	*pathdir, *base_name;
        int     n_samps, Ysz, Xsz, n_episodes, min_pass_duration;
	EPISODE	*episodes;
{
        char    *strcat(),*strcpy();
        FILE    *fopen();
        void    find_paths();
	void	calc_z();
	void	calc_locomotion();

        FILE    *tfp,*zfp;
        int     i,j,s, cent_y, cent_x, start_s, stop_s;
        char	out_file[256],temp_name[256];
        INT_MAP hotspot;
        unsigned int    x, y, x0,y0,x1,y1;
        int     good_samps, pass_duration;
	int	obs;
	float	speed, linearity, path_length, exp;
	double	z;


        (void)sprintf(out_file,"%s/%s.T",pathdir, base_name);
        if((tfp = fopen(out_file,"w")) == NULL){
                (void)fprintf(stderr,"Can't open output file %s\n",out_file);
                exit(1);
        }
	(void)sprintf(temp_name,"%s/%s.Tz",pathdir,base_name);
        if((zfp = fopen(temp_name,"w")) == NULL){
                (void)fprintf(stderr,"Can't open output file %s\n",temp_name);
                exit(1);
        }

	for(i=0; i < n_episodes; i++){
                if(episodes[i].stop > n_samps)
                        continue;
                start_s = episodes[i].start;
                stop_s = episodes[i].stop;
                pass_duration = stop_s - start_s;
                if(pass_duration < min_pass_duration){
                	(void)fprintf(zfp,"%d\t%8.5f\t%8.5lf\t%8.5f\t%d\t%d\t%8.3f\t%8.3f\t%8.3f\n",0, 0.0, 0.0, 0.0, 0, 0, 0.0,0.0, 0.0);
                	continue;
		}
                if(check_sampling(pos, start_s, stop_s, &good_samps, time, min_samples)){
                	(void)fprintf(zfp,"%d\t%8.5f\t%8.5lf\t%8.5f\t%d\t%d\t%8.3f\t%8.3f\t%8.3f\n",0, 0.0, 0.0, 0.0, 0, 0, 0.0,0.0, 0.0);
                        continue;
		}

		// find the first and last good samples in the episode for the trajectory header
		for(s=start_s; s <= stop_s; s++){
			y0 = (int)(pos[s].y);
			x0 = (int)(pos[s].x);
                        if(y0 && x0)
				break;
		}
		for(s=stop_s; s >= start_s; s--){
			y1 = (int)(pos[s].y);
			x1 = (int)(pos[s].x);
                        if(y1 && x1)
				break;
		}
		// write the trajectory header
                fprintf(tfp,"%d\t%d  %d\t%d  %d\t%d  %d\n",good_samps,y0,x0,0,0,y1,x1);
			
		// write the trajectory
		for(s=start_s; s < stop_s; s++){
			y = (int)(pos[s].y);
			x = (int)(pos[s].x);
                        if(y && x){
               		         (void)fprintf(tfp,"%d\t%d\t%d\t%2.6f\n",y,x,spk[s],rate[y][x]/ (double)samps_per_sec);
                         }
                }

		// calculate the measures of the trajectory
		calc_z(start_s, pass_duration, pos, spk, rate, &obs, &exp, &z, samps_per_sec);
		// calculate the measures of the trajectory
		calc_locomotion(i, pass_duration, pos, samps_per_sec, &path_length, &speed, &linearity);

		// write the measures to the z file
                (void)fprintf(zfp,"%d\t%8.5f\t%8.5lf\t%8.5f\t%d\t%d\t%8.3f\t%8.3f\t%8.3f\n",obs, exp, z, grand_rate, start_s, pass_duration, path_length, speed, linearity);

	}		

        fclose(tfp);
        fclose(zfp);
        return;
}
