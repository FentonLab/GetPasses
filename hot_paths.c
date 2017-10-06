/*
	outfile is infile.f#T , where f# is field number. The file is :
	n_samples yin xin yhot xhot yout xout
	y x s rt ... for each path. rt is the session rt for the pixel. The path begins and ends in the field. The output is the sequence of positions and spikes and the expected rate (y x s rt) that include a pass through the "hot_spot" (field center)*/

#include		<stdio.h>
#include		<stdlib.h>
#include		<time.h>
#include		<math.h>
#include		"apparat.h"
#include		"cell.h"
#include		"translate.h"
#include		"printer.h"
#include		"hot_paths.h"

#define	TWOPI		6.283185307
#define			ARM_DIST	6.0

SESS		sess;
DETECTS		dets;
CELL		cell;
APPAR		app;
MY_ERR		my_err;

void		do_err();
void		clr_err();

int		reduce_flag = FALSE;	/* if set, reduce position resolution 2:1 */
int		two_reducts = FALSE;	/* reduce postion resolution 4:1; d and r flags cannot both be
					   set */
int		which_cell  = 1;	/* do cell 1 if no other is indicated */
int		tm_out_flag = FALSE;	/* means time out bit will be ignored */
int		gray_flag   = FALSE;	/* print peak rate pixels for fields in gray */
int		eight_flag  = FALSE;	/* Default - use 6 cols for sample map, otherwise 8 cols */
int		tek_flag    = FALSE;	/* Print to tek except if junk_flag  */
int		alpha_flag  = FALSE;
int		pg_1_maps   = FALSE;
int		paper_sz    = SMALL;
int		French_flag = FALSE;
str256		temp_nrs;
str256		prog_name;
str256		menu_dir;
str256		my_line;
str256		command,change_dir,data_dir;
float		samp_per_sec = 60.0;
int		samp_freq = 60;
int		which_spot = DAT_DIR;
int             spot_sz = 1;
int             path_length = 300;
int             min_path_length = 60;
int		pedit_file = 0;
int		start_samp = 0;
int		mid_times = 0;
int		show_times = 0;
int		input_paths = 0;
int		arbitrary_paths = 0;
int		arbitrary_length = 0;
int		min_dist = 0;
int		samps_after = 0;
int		samps_before = 0;
int		eight_arm = 0;

extern	int	optind;
extern	char	*optarg;
extern	double	atof();

main(argc,argv)
int	argc;
char	*argv[];
{
	int	c;
	char	*strcat(),*strcpy(),*getenv();
	void	instruct();
	void	do_basic_stuff_2B();
	void	exit();
	FILE	*fopen(),*ifp;
	void	find_hot_paths();
	void	write_posfile();
	void	*calloc();
	int	test_run_speed();
	int	check_sampling();

	int	part_sess = 0, min_time = 15;
	str256	base_name,input_times_file,pathdir,miscdir;
	int	i,j, k, temp,n_paths = 0, max_samps;
	int	integ_time = 0, good_samps;
	PATH	*pre_def_path;
	char	*samps_taken;

        struct timeval  tp;
        struct timezone tzp;
        USI     seed[3];
 	int     gettimeofday();
        double  drand48();
        USI     *seed48();

	strcpy(prog_name,argv[0]);
	strcpy(temp_nrs,"/users/angle/tmp/temp_nrs");

        if(argc == 1)
                instruct();

        cell.fld_thr = EPS;
        cell.cell_nr = which_cell;
        cell.unit_nr = 0;

	arbitrary_length = 3 * samp_freq;
        while((c = getopt(argc,argv,"a:b:c:d:ef:i:j:k:l:mn:o:prs:t:w")) != EOF)
                switch(c){
			case	'a':	samps_after = atoi(optarg);
					break;
			case	'b':	samps_before = atoi(optarg);
					break;
			case	'c':	min_time = atoi(optarg);
					break;
			case	'd':	part_sess = atoi(optarg);
					part_sess *= 3600;
					break;
			case	'e':	eight_arm = 1;
					break;
			case	'f':	samp_freq = atoi(optarg);
					samp_per_sec = (float)samp_freq;
					break;
			case	'i':	input_paths= 1;
					strcpy(input_times_file,optarg);
					break;
			case	'j':	arbitrary_paths = atoi(optarg);
					break;
			case	'k':	arbitrary_length = atoi(optarg);
					break;
			case	'l':	min_dist = atoi(optarg);
					break;
			case	'm':	mid_times = TRUE;
					break;
                        case    'n':    min_path_length = atoi(optarg);
                                        break;
			case	'o':	start_samp = atoi(optarg);
					break;
			case	'p':	pedit_file = 1;
					break;
                        case    'r':    reduce_flag = TRUE;
                                        break;
                        case    's':    spot_sz = atoi(optarg);
                                        break;
                        case    't':    cell.fld_thr = atof(optarg);
                                        break;
                        case    'w':   	show_times = TRUE;
                                        break;

                        default:        clr_err();
                                        sprintf(my_err,"Unknown option on command line.",argv[0]);
                                        do_err("main",1);
                }

        if(argc - optind != 1){
                clr_err();
                sprintf(my_err,"Call with options and a one-spot (2 bytes / samp) data file in the DATA_DIR.");
                do_err("main",3);
        }
	if(input_paths && arbitrary_paths)
		instruct();

	if(*(strcpy(pathdir,getenv("PATH_DIR"))) == NULL){
		(void)fprintf(stderr,"PATH_DIR not defined in environment.\n");
		exit(2);
	}
	if(*(strcpy(miscdir,getenv("MISC_DIR"))) == NULL){
		(void)fprintf(stderr,"MISC_DIR not defined in environment.\n");
		exit(2);
	}

	integ_time = (int)(0.25 * samp_per_sec);
	(void) strcpy(sess.file_name,argv[optind]);

	do_basic_stuff_2B(reduce_flag,tm_out_flag,part_sess);
	do_rate_stuff();
	cell.n_fields = get_field_map(cell.al_rt,cell.fld_map,cell.fld_thr);
	cell.n_fields = prune_field_list(cell.fld_map,cell.n_fields,MIN_FLD_SIZ);
	do_field_stuff();
	set_pr_fields();
	do_rate_map();
	make_time_map(cell.al_tm,cell.sample_time);

	if(input_paths){
		(void)strcpy(base_name,miscdir);
        	(void)strcat(base_name,input_times_file);
        	if((ifp = fopen(base_name,"r")) == NULL){
                	(void)fprintf(stderr,"Can't open path defining samples file %s\n",base_name);
                	exit(1);
		}
		n_paths = 0;
		while((fscanf(ifp,"%d%d",&temp,&temp)!= EOF))
			n_paths++;
		if((pre_def_path = (PATH *)calloc((size_t)n_paths,sizeof(PATH))) == NULL){
			fprintf(stderr,"Can't calloc pre_def_paths\n");
			return;
		}
		rewind(ifp);
		for(i=0; i < n_paths;i++)
			fscanf(ifp,"%d%d",&(pre_def_path[i].start), &(pre_def_path[i].stop));
        	fclose(ifp);
	}

	if(arbitrary_paths){
        	(void)gettimeofday(&tp,&tzp);
        	seed[0] = (USI)tp.tv_usec % MAX_SHORT;
        	(void)seed48(seed);

		n_paths = arbitrary_paths;
		if((pre_def_path = (PATH *)calloc((size_t)n_paths,sizeof(PATH))) == NULL){
			fprintf(stderr,"Can't calloc pre_def_paths\n");
			return;
		}
		max_samps = sess.file_size / BYTES_PER_SAMP_2B; 
		if((samps_taken = (char *)calloc((size_t)max_samps,sizeof(char))) == NULL){
			fprintf(stderr,"Can't calloc samps_taken\n");
			return;
		}
		for(k=0; k < max_samps; k++)
			samps_taken[k] = 0;
		max_samps -= arbitrary_length - 1;
		for(i=0; i < n_paths; i++){
			while(1){
				TRY_AGAIN:;
				pre_def_path[i].start = (int)(drand48() * (double)max_samps);
				for(j=0, k=pre_def_path[i].start; j < arbitrary_length; j++, k++)
					if(samps_taken[k])
						goto TRY_AGAIN;
					
				if(check_sampling(pre_def_path[i].start,pre_def_path[i].start + arbitrary_length, &good_samps, cell.al_tm, min_time))
					goto TRY_AGAIN;
				if(good_samps < arbitrary_length)
					goto TRY_AGAIN;

				if(min_dist){		/* test for min running speed */
					if(!test_run_speed(pre_def_path[i].start, pre_def_path[i].start + arbitrary_length, integ_time, (double)min_dist))
						goto TRY_AGAIN;
				}

				for(j=0, k=pre_def_path[i].start; j < arbitrary_length; j++, k++)
					samps_taken[k] = 1;
				pre_def_path[i].stop = k;
				break;
			}
		}
	}

	if(pedit_file)
		write_posfile(miscdir);
	find_hot_paths(pathdir,miscdir,pre_def_path,n_paths,min_time);
	return;
}

void instruct()
{
        void    exit();

        fprintf(stderr,"%s: Call with options and a 2 byte/samp KM data file\n",prog_name);
	fprintf(stderr,"Use to select and/or define paths\n");
	fprintf(stderr,"INPUT: KM file in DATA_DIR\n");
	fprintf(stderr,"OUTPUT: 'trajectory' file (.T) in PATH_DIR\n");
	fprintf(stderr,"EXPECTS a MISC_DIR\n");
        fprintf(stderr,"Options:\n");
	fprintf(stderr,"\t-a: Define path # samples after the field entry.\n");
	fprintf(stderr,"\t-b: Define path # samples before the field entry.\n");
	fprintf(stderr,"\t-c: Min number of samples that all positions in a path must be visited. Default 15.\n");
	fprintf(stderr,"\t-d: Set session duration in minutes (default: all samples).\n");
	fprintf(stderr,"\t-e: TEMPORARY write to stdout arm (1-8) y x latency (samps) to go to an arm of the 8-arm maze before entering and after leaving the field\n");
	fprintf(stderr,"\t-f: Set sampling frequency (default: 60 samps per sec).\n");
	fprintf(stderr,"\t-i: Use a file of 60 Hz sample pairs to select  paths. Give file name. Expected in MISC_DIR.\n");
	fprintf(stderr,"\t-j: Define N randomly selected paths of length specified by option 'k'.\n");
	fprintf(stderr,"\t-k: Set path length for randomly selected paths. (Default = 3 x sampling freq)\n");
	fprintf(stderr,"\t-l: if -k: only use paths when the rat is running at a min dist / sec.\n");
        fprintf(stderr,"\t-m: Output to MISC_DIR a file (infilem) of mid-path times (60 Hz sample times)\n");
        fprintf(stderr,"\t-n: set the mininum path length (in 60 Hz sample times). default = 60\n");
	fprintf(stderr,"\t-o: Not implemented yet. Offset the start sample. To put in register with a movie\n");
        fprintf(stderr,"\t-p: Output to MISC_DIR:\n");
	fprintf(stderr,"\t    the list of positions (infilePOS) values are scaled for 64X64 and optionally compressed\n");
	fprintf(stderr,"\t    the edit list (infileE) in 30 Hz frame units.\n");
	fprintf(stderr,"\t    the spikes file (infileS) in 30 Hz time stamp units.\n");
	fprintf(stderr,"\t    the field file (infile.#FLD). Hot spot indicated by adding 1000 to pixel values.\n");
        fprintf(stderr,"\t-r:  set reduce flag (compress 2:1)\n");
        fprintf(stderr,"\t-s: Hot spot size (radius in pixels)\n");
        fprintf(stderr,"\t-t: Set threshold for boundary between zero and finite rate\n");
        fprintf(stderr,"\t-w: Output to MISC_DIR a file (infilew) of path start and stop (60 Hz) sample times.\n");
        exit(0);
}

char	*zap_newline();

#define		PR_1			fprintf(printer,"%s\n",my_line)
#define		N_WIDE_IN_LN		80
#define		LINES_PER_PG		60

int	pr_fields;
str256	pr_fld_warn;

set_pr_fields()
{
	if(cell.n_fields > MAX_PRINT_FIELDS){
		pr_fields = MAX_PRINT_FIELDS;
		sprintf(pr_fld_warn,"Data are shown only for the %d largest fields",MAX_PRINT_FIELDS);
		}
	else{
		pr_fields = cell.n_fields;
		pr_fld_warn[0] = '\000';
		}
}

char	*itoa(n)
int	n;
{
	static	char	val[20];

	sprintf(val,"%d",n);

	return val;
}

void    write_posfile(miscdir)
	str256	miscdir;
{
	FILE	*fopen();
        char    *strcat(),*strcpy();

        FILE    *pfp,*sfp;
        int     i,b,samp;
        str256  out_file;
        uint    byte[BYTES_PER_SAMP_2B];
        uint    x,y,spks;
        
        (void)sprintf(out_file,"%s%sPOS",miscdir,sess.file_name);
        if((pfp = fopen(out_file,"w")) == NULL){
                (void)fprintf(stderr,"Can't open positions file %s\n",out_file);
                exit(1);
        }
        (void)sprintf(out_file,"%s%sS",miscdir,sess.file_name);
        if((sfp = fopen(out_file,"w")) == NULL){
                (void)fprintf(stderr,"Can't open spikes file: %s\n",out_file);
                exit(1);
        }
        
	spks = 0;
	for(i=0,samp = 0;i < sess.file_size;i += BYTES_PER_SAMP_2B,samp++){
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
		fprintf(pfp,"%d\t%d\t%d\n",samp,y,x);
                
		spks += byte[1] >> 6;
                spks += (byte[0] >> 6) << 2;
		if(spks)
			if(samp % SAMPS_PER_FRAME){
				fprintf(sfp,"%d\t%d\n",samp / SAMPS_PER_FRAME,spks);
				spks = 0;
			}
	}
	fclose(pfp);
	fclose(sfp);
	return;
}

void    find_hot_paths(pathdir, miscdir,pre_def_path,n_paths,min_time)
	str256	pathdir, miscdir;
	PATH	*pre_def_path;
	int	n_paths,min_time;
{
        FILE    *ofp,*efp,*ffp,*tmfp,*twfp;
        int     i,j,fld_nr;
        str256  out_file,base_name,temp_name;
        u_map   hotspot;
        int     x_start,y_start,x_stop,y_stop;

        void    zero_usi();
        char    *strcat(),*strcpy();
	FILE	*fopen();
	void	use_path_defs(), find_paths();


        (void)strcpy(base_name,miscdir);
        (void)strcat(base_name,sess.file_name);

        for(fld_nr=0;fld_nr < cell.n_fields;fld_nr++){
		if(pedit_file){
        		(void)sprintf(out_file,"%s.%dE",base_name,fld_nr);
        		if((efp = fopen(out_file,"w")) == NULL){
               		 	(void)fprintf(stderr,"Can't open edit file %s\n",out_file);
               		 	exit(1);
			}
        	}
                (void)sprintf(out_file,"%s%s.%dT",pathdir,sess.file_name,fld_nr);
                if((ofp = fopen(out_file,"w")) == NULL){
                        (void)fprintf(stderr,"Can't open output file %s\n",out_file);
                        exit(1);
                }
                if(mid_times){
                	(void)sprintf(temp_name,"%s.%dTm",base_name,fld_nr);
			if((tmfp = fopen(temp_name,"w")) == NULL){
                        	(void)fprintf(stderr,"Can't open output file %s\n",temp_name);
                        	exit(1);
			}
                }
                
                if(show_times){
                	(void)sprintf(temp_name,"%s.%dTw",base_name,fld_nr);
			if((twfp = fopen(temp_name,"w")) == NULL){
                        	(void)fprintf(stderr,"Can't open output file %s\n",temp_name);
                        	exit(1);
			}
                }
		zero_usi(hotspot);
                y_start = cell.fields[fld_nr].y_ctr - spot_sz + 1;
                y_stop = cell.fields[fld_nr].y_ctr + spot_sz - 1;
                x_start = cell.fields[fld_nr].x_ctr - spot_sz + 1;
                x_stop = cell.fields[fld_nr].x_ctr + spot_sz - 1;
		
		for(i= y_start;i <= y_stop;i++)
                        for(j= x_start;j <= x_stop;j++)
                                hotspot[i][j] = (cell.fld_map[i][j] == fld_nr+1);
		if(pedit_file){
        		(void)sprintf(out_file,"%s.%dFLD",base_name,fld_nr);
        		if((ffp = fopen(out_file,"w")) == NULL){
               		 	(void)fprintf(stderr,"Can't open field file %s\n",out_file);
               		 	exit(1);
			}
			for(i= 0;i < YSZ;i++)
                       		for(j= 0;j < XSZ;j++){
					if((i >= y_start) && (i <= y_stop) && (j >= x_start) && (j <= x_stop)){
						fprintf(ffp,"%d\t%d\n",i+HOT_SPOT,j+HOT_SPOT);
						continue;
					}
                                	if(cell.fld_map[i][j] == fld_nr+1)
						fprintf(ffp,"%d\t%d\n",i,j);
				}

			fclose(ffp);
		}
		if(n_paths)
			use_path_defs(ofp,efp,tmfp,twfp,pre_def_path,n_paths);
		else
			find_paths(fld_nr,hotspot,ofp,efp,tmfp,twfp,min_time);

        	
		fclose(ofp);
		if(pedit_file)
			fclose(efp);
	}
        return;
}

int	enter_field(c,i,y0,x0,y1,x1)
	int	c,*i;
	uint	*y0,*x0,*y1,*x1;
	{
	uint	byte[BYTES_PER_SAMP_2B];
	uint	 s,temp,y=0,x=0,b;
	int	n_samps = 0;
	
	s = temp = *i;

	do{
		*i = s;	
		*y0 = y;
		*x0 = x;
		n_samps++;
		s-= BYTES_PER_SAMP_2B;
        	for(b=0;b<BYTES_PER_SAMP_2B;b++)
                	byte[b] = sess.data[s + b];
		if(reduce_flag){
                       	y = (byte[0] & 63) >> 1;
                       	x = (byte[1] & 63) >> 1;
			if(y && x){
				y += 16;
				x += 16;
			}
		}else{
                       	y = (byte[0] & 63);
                        x = (byte[1] & 63);
		}
	}while((cell.fld_map[y][x] == c+1) || !(y && x));
	s = temp;
	do{
		*y1 = y;
		*x1 = x;
		s+= BYTES_PER_SAMP_2B;
		n_samps++;
        	for(b=0;b<BYTES_PER_SAMP_2B;b++)
                	byte[b] = sess.data[s + b];
		if(reduce_flag){
                       	y = (byte[0] & 63) >> 1;
                       	x = (byte[1] & 63) >> 1;
                       if(y && x){
                                y += 16;
                                x += 16;
                        }
		}else{
                       	y = (byte[0] & 63);
                        x = (byte[1] & 63);
		}
	}while((cell.fld_map[y][x] == c+1) || !(y && x));
	if((*i -= samps_before) < 0)
		*i = 0;
	return(n_samps + samps_before + samps_after -1);
}
			
void	find_paths(fld_nr,hotspot,ofp,efp,tmfp,twfp,min_time)
	int	fld_nr,min_time;
        u_map   hotspot;
	FILE	*ofp,*efp,*tmfp,*twfp;

{
        int     i,b,s;
        uint    x,y;
        uint    x0,y0,x1,y1;
        uint    tx,ty;
	int	y8, x8;
	int	spks;
        uint    byte[BYTES_PER_SAMP_2B];
	int	in_samp,out_samp, good_samps, d_samps, bin, i8;
	double	angle, hypot(), atan2(), d8;


        for(i=0;i < sess.file_size;i += BYTES_PER_SAMP_2B){
        	for(b=0;b<BYTES_PER_SAMP_2B;b++)
        		byte[b] = sess.data[i + b];

		if(reduce_flag){
       		        y = (byte[0] & 63) >> 1;
       		        x = (byte[1] & 63) >> 1;
       		        if(y && x){
       	        		y += 16;
       	                	x += 16;
       	        	}
		}else{
       	        	y = (byte[0] & 63);
       	        	x = (byte[1] & 63);
		}
	
       	 	if(hotspot[y][x]){
			path_length = enter_field(fld_nr,&i,&y0,&x0,&y1,&x1);
			in_samp = i / BYTES_PER_SAMP_2B;
			out_samp = in_samp + path_length - 1;
			if((i + (path_length * BYTES_PER_SAMP_2B)) > sess.file_size)
				break;
			if(path_length < min_path_length)
				break;
			if(check_sampling(in_samp, out_samp, &good_samps, cell.al_tm, min_time))
				break;	

			if(mid_times)
				fprintf(tmfp,"%d\n",in_samp + ((out_samp - in_samp) / 2));  

			if(show_times)
				fprintf(twfp,"%d\t%d\n",in_samp,out_samp);
			
			(void)fprintf(ofp,"%d\t%d  %d\t%d  %d\t%d  %d\n",good_samps,y0,x0,y,x,y1,x1);

			i8 = i;
                	for(s=0;s < path_length;s++,i+= BYTES_PER_SAMP_2B){
                        	for(b=0;b<BYTES_PER_SAMP_2B;b++)
                			byte[b] = sess.data[i + b];
				if(reduce_flag){
                     			ty = ((byte[0] & 63) >> 1);
                      			tx = ((byte[1] & 63) >> 1);
                     			if(ty && tx){
                             			ty += 16;
                              			tx += 16;
                      			}
				}else{
                     			ty = (byte[0] & 63);
                      			tx = (byte[1] & 63);
				}
				if(ty && tx){
                        		spks = byte[1] >> 6;
                        		spks += (byte[0] >> 6) << 2;
                        		(void)fprintf(ofp,"%d\t%d\t%d\t%2.6f\n",ty,tx,spks,cell.al_rt[ty][tx]/ (double)samp_per_sec);
				}
               		}
			if(pedit_file)
				fprintf(efp,"%d\t%d\n",in_samp / SAMPS_PER_FRAME,out_samp / SAMPS_PER_FRAME);
			if(eight_arm){
				d_samps = 0;
				while (i8 > 0){
                        		for(b=0;b<BYTES_PER_SAMP_2B;b++)
                				byte[b] = sess.data[i8 + b];
					d_samps++;
					if(reduce_flag){
               		      			y8 = (int)((byte[0] & 63) >> 1);
               		       			x8 = (int)((byte[1] & 63) >> 1);
               		      			if(y8 && x8){
               		              			y8 += 16;
               		               			x8 += 16;
               		       			}
					}else{
               		      			y8 = (int)(byte[0] & 63);
               		       			x8 = (int)(byte[1] & 63);
					}
					if(y8 && x8){
						d8 = hypot((double)((int)y0 - y8), (double)((int)x0 - x8));
						if(d8 > ARM_DIST){
							angle = atan2((double)(32 - y8), (double)(x8 - 32));
							angle -= (TWOPI / 16.0);
                                			if(angle < 0.0)
                                        			angle += TWOPI;
                                			if(angle > TWOPI)
                                        			angle -= TWOPI;
                                			angle /= TWOPI;
                                			angle *= 8.0;
                                			bin =  1 + (int)(angle + 0.5);
							if(bin == 9)
								bin = 1;
							printf("%d\t%d\t%d\t%d\t\t",bin, y8, x8, d_samps);				
							break;
						}
					}
					i8 -= BYTES_PER_SAMP_2B;
				}
				if(d8 <= ARM_DIST)
					printf("%d\t%d\t%d\t%d\t\t",9, 0, 0, 0);				
				i8 = i;
				d_samps = 0;
				while (i8  < sess.file_size){
                        		for(b=0;b<BYTES_PER_SAMP_2B;b++)
                				byte[b] = sess.data[i8 + b];
					d_samps++;
					if(reduce_flag){
               		      			y8 = (int)((byte[0] & 63) >> 1);
               		       			x8 = (int)((byte[1] & 63) >> 1);
               		      			if(y8 && x8){
               		              			y8 += 16;
               		               			x8 += 16;
               		       			}
					}else{
               		      			y8 = (byte[0] & 63);
               		       			x8 = (byte[1] & 63);
					}
					if(y8 && x8){
						d8 = hypot((double)((int)y1 - y8), (double)((int)x1 - x8));
						if (d8 > ARM_DIST){
							angle = atan2((double)(32 - y8), (double)(x8 - 32));
							angle -= (TWOPI / 16.0);
                                			if(angle < 0.0)
                                        			angle += TWOPI;
                                			if(angle > TWOPI)
                                        			angle -= TWOPI;
                                			angle /= TWOPI;
                                			angle *= 8.0;
                                			bin = 1 + (int)(angle + 0.5);
							if(bin == 9)
								bin = 1;
							printf("%d\t%d\t%d\t%d\n",bin, y8, x8, d_samps);				
							break;
						}
					}
					i8 += 2;
				}
				if(d8 <= ARM_DIST)
					printf("%d\t%d\t%d\t%d\n",9, 0, 0, 0);				
			}
				
		}
	}
	return;
}
void	use_path_defs(ofp,efp,tmfp,twfp,pre_def_path,n_paths)
	FILE	*ofp,*efp,*tmfp,*twfp;
	PATH	*pre_def_path;
	int	n_paths;
{
        int     i,b,s,p;
        uint    x0,y0,x1,y1;
        uint    tx,ty;
	int	spks;
        uint    byte[BYTES_PER_SAMP_2B];
	int	n_samps;
	int	start_s,stop_s,start_b,stop_b;

	n_samps = sess.file_size / BYTES_PER_SAMP_2B;

	for(p=0;p < n_paths;p++){
		if(pre_def_path[p].stop > n_samps)
			continue;
		start_s = pre_def_path[p].start;
		stop_s = pre_def_path[p].stop;
		start_b = start_s * BYTES_PER_SAMP_2B;
		stop_b = stop_s * BYTES_PER_SAMP_2B;

		path_length = stop_s - start_s + 1;

		i = start_b;
                for(b=0;b<BYTES_PER_SAMP_2B;b++)
                	byte[b] = sess.data[i + b];
		if(reduce_flag){
                	y0 = ((byte[0] & 63) >> 1);
                     	x0 = ((byte[1] & 63) >> 1);
                    	if(y0 && x0){
                           	y0 += 16;
                            	x0 += 16;
                     	}
		}else{
                    	y0 = (byte[0] & 63);
                     	x0 = (byte[1] & 63);
		}
		i = stop_b;
                for(b=0;b<BYTES_PER_SAMP_2B;b++)
                	byte[b] = sess.data[i + b];
		if(reduce_flag){
                	y1 = ((byte[0] & 63) >> 1);
                     	x1 = ((byte[1] & 63) >> 1);
                    	if(y1 && x1){
                           	y1 += 16;
                            	x1 += 16;
                     	}
		}else{
                    	y1 = (byte[0] & 63);
                     	x1 = (byte[1] & 63);
		}

		if(mid_times)
			fprintf(tmfp,"%d\n",start_s + ((stop_s - start_s) / 2));  

		if(show_times)
			fprintf(twfp,"%d\t%d\n",start_s,stop_s);
			
		(void)fprintf(ofp,"%d\t%d  %d\t%d  %d\t%d  %d\n",path_length,y0,x0,0,0,y1,x1);

               	for(i = start_b, s=0;s < path_length;s++,i+= BYTES_PER_SAMP_2B){
                        for(b=0;b<BYTES_PER_SAMP_2B;b++)
                		byte[b] = sess.data[i + b];
			if(reduce_flag){
                     		ty = ((byte[0] & 63) >> 1);
                      		tx = ((byte[1] & 63) >> 1);
                     		if(ty && tx){
                             		ty += 16;
                              		tx += 16;
                      		}
			}else{
                     		ty = (byte[0] & 63);
                      		tx = (byte[1] & 63);
			}
                        spks = byte[1] >> 6;
                        spks += (byte[0] >> 6) << 2;
                        (void)fprintf(ofp,"%d\t%d\t%d\t%2.2f\n",ty,tx,spks,cell.al_rt[ty][tx]/ (double)samp_per_sec);
               		
			if(pedit_file)
				fprintf(efp,"%d\t%d\n",start_s / SAMPS_PER_FRAME,stop_s / SAMPS_PER_FRAME);
		}
	}
	return;
}
int	check_sampling(in_samp, out_samp, good_samps,time, min_time)
	int	in_samp, out_samp, *good_samps, min_time;
	TIME	time;
{
	int	i, b, y, x, start, stop;
        uint    byte[BYTES_PER_SAMP_2B];
	
	*good_samps = 0;
	start = in_samp * BYTES_PER_SAMP_2B;
	stop = out_samp * BYTES_PER_SAMP_2B;

        for(i= start;i <= stop;i += BYTES_PER_SAMP_2B){
		for(b=0;b <BYTES_PER_SAMP_2B;b++)
               		byte[b] = sess.data[i + b];

               	if(reduce_flag){
                 	y = (byte[0] & 63) >> 1;
                        x = (byte[1] & 63) >> 1;
                        if(y && x){
                               	y += 16;
                               	x += 16;
                        }
                }else{
                        y = (byte[0] & 63);
                        x = (byte[1] & 63);
                }
		if(y && x){
			(*good_samps)++;
			if(time[y][x] < min_time)
				return(1);
		}
	}	
	return 0;
}
	
int	test_run_speed(start, stop, integ_time, min_distance)
	int	start, stop, integ_time;
	double	min_distance;
{
	int	i, b, y, x, y0 = 0, x0 = 0;
	double hypot(), dist;
        uint    byte[BYTES_PER_SAMP_2B];

	min_dist /= 4.0;
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
