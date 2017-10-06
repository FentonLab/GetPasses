#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  "tsf.h"


static neighbors_1();
static int	compare_1();

mk_field_map(rate, time, fm, thr, Ysz, Xsz)
FLOAT_MAP       rate;
INT_MAP         time,fm;
float           thr;
int             Ysz, Xsz;

{
	int	i, j, fld_num = 1;

        for(i=0; i < Ysz; i++)
                for(j=0; j < Xsz; j++)
                        if(time[i][j] != -1)
                                fm[i][j] = 0;
                        else
                                fm[i][j] = -1;

	for (i = 1; i < Ysz - 1; i++)
		for (j = 1; j < Xsz - 1; j++)
			if (rate[i][j] >= thr)
				if (fm[i][j] == 0) {
					neighbors_1(i, j, rate, fm, thr, fld_num);
					fld_num++;
				}
	return(fld_num - 1);
}


neighbors_1(i, j, r, fm, thr, fld_num)
int	i, j, fld_num;
FLOAT_MAP	r;
float	thr;
INT_MAP	fm;

{
	int	itmp, jtmp;

	fm[i][j] = fld_num;
	itmp = i;

	for (i -= 1; i < itmp + 2; i += 2){
		if(i < 0)
			continue;
		if (fm[i][j] == 0)
			if (r[i][j] >= thr)
				neighbors_1(i, j, r, fm, thr, fld_num);
	}

	i = itmp;
	jtmp = j;

	for (j -= 1; j < jtmp + 2; j += 2){
		if(j < 0)
			continue;
		if (fm[i][j] == 0)
			if (r[i][j] >= thr)
				neighbors_1(i, j, r, fm, thr, fld_num);
	}
}

int	prune_field_list(fmap, n_flds, min_siz,Ysz,Xsz)
INT_MAP	fmap;
int	n_flds, min_siz,Ysz,Xsz;
{
	void	 *calloc();
	int	i, j, n;
	int	tmp, good_flds;
	INT_MAP	t_map;
	SORT_FLDS	*to_sort;

	to_sort = (SORT_FLDS *) calloc((unsigned)n_flds, sizeof(SORT_FLDS));

	for (n = 1; n <= n_flds; n++)
		to_sort[n-1].f_nr = n;

	for (i = 0; i < Ysz; i++)
		for (j = 0; j < Xsz; j++)
			if (fmap[i][j] > 0)
				to_sort[fmap[i][j] - 1].f_sz++;

	qsort((char * ) to_sort, (unsigned)n_flds, sizeof(SORT_FLDS), compare_1);

	n = 0;
	while (n < n_flds && to_sort[n].f_sz >= min_siz) 
		n++;

	good_flds = n;

	/* This next loop need not be done if n == n_flds, but what the heck */

	for (i = 0; i < Ysz; i++)
		for (j = 0; j < Xsz; j++)
			if ((tmp = fmap[i][j]) > 0)
				for (n = good_flds; n < n_flds; n++)
					if (tmp == to_sort[n].f_nr) {
						fmap[i][j] = 0;
						break;
					}

	zero_int(t_map,Ysz,Xsz);

	for (i = 0; i < Ysz; i++)
		for (j = 0; j < Xsz; j++) {
			for (n = 0; n < good_flds; n++)
				if (fmap[i][j] == to_sort[n].f_nr) {
					t_map[i][j] = n + 1;
					break;
				}
		}

	copy_int(fmap, t_map,Ysz,Xsz);
	return(good_flds);
}


static compare_1(a, b)
SORT_FLDS	*a, *b;
{
	if (a->f_sz > b->f_sz) 
		return(-1);
	if (a->f_sz < b->f_sz) 
		return(1);
	return(0);
}


void	do_field_stuff(fm, Ysz, Xsz, cent_y, cent_x)
FRAME_MAP	*fm;
int	Ysz, Xsz;
float	cent_y, cent_x;
{
	void	 * calloc();
	int	f, fnr, i, j;
	void	cent_gravity(), get_pol_center();

	fm->oof_time = 0;
	fm->oof_spks = 0;
	fm->if_spks = 0;
	fm->if_time = 0;
	for (i = 0; i < Ysz; i++){
		for (j = 0; j < Xsz; j++){
			if(fm->field_map[i][j] == 0){
				if(fm->time_map[i][j] > 0){
					fm->oof_time += fm->time_map[i][j];
					fm->oof_spks += fm->spk_map[i][j];
				}
			}else{
				 if (fm->field_map[i][j] > 0){
					fm->if_time += fm->time_map[i][j];
					fm->if_spks += fm->spk_map[i][j];
				}
			}
		}
	}
	fm->oof_time /= fm->samps_per_sec;
	fm->if_time /= fm->samps_per_sec;

	if(fm->oof_time)
		fm->oof_rate = (float)fm->oof_spks / (float)fm->oof_time;

	for (f = 0, fnr = 1; f < fm->n_fields; f++, fnr++) {
		fm->field[f].size = 0;
		fm->field[f].spks = 0;
		fm->field[f].time = 0;
		fm->field[f].oof_act_pix = 0;
		fm->field[f].oof_spks = 0;
		fm->field[f].oof_time = 0;
		for (i = 0; i < Ysz; i++){
			for (j = 0; j < Xsz; j++){
				if (fnr == fm->field_map[i][j]) {
					fm->field[f].size++;
					fm->field[f].time += fm->time_map[i][j];
					fm->field[f].spks += fm->spk_map[i][j];
				}else {
					if(fm->time_map[i][j] > 0){
						fm->field[f].oof_time += fm->time_map[i][j];
						if(fm->spk_map[i][j] > 0){
							fm->field[f].oof_act_pix++;
							fm->field[f].oof_spks += fm->spk_map[i][j];
						}
					}
				}
			}
		}
		fm->field[f].time /= fm->samps_per_sec;	
		fm->field[f].oof_time /= fm->samps_per_sec;	

		fm->field[f].num = fnr;
		fm->field[f].rate = (float) fm->field[f].spks / (float)fm->field[f].time;
		fm->field[f].rate *= (float) PRAHA_SAMPS_PER_SEC;
		if(fm->field[f].oof_time > 0)
			fm->field[f].oof_rate = (float)fm->field[f].oof_spks / (float)fm->field[f].oof_time;
		else
			fm->field[f].oof_rate = -1.0;
		if(fm->field[f].oof_rate > 0.0){
			fm->field[f].s2n = fm->field[f].rate / fm->field[f].oof_rate;
			fm->field[f].concent = (fm->field[f].rate / fm->field[f].size);
		}else{
			fm->field[f].s2n = -1.0;
			fm->field[f].concent = -1.0;
		}

		get_field_center(&(fm->field[f]), fnr, fm->field_map, fm->time_map, fm->spk_map,Ysz, Xsz, cent_y, cent_x);
		cent_gravity(fm->rate_map,fm->field_map,fnr,&(fm->field[f].cg_y),&(fm->field[f].cg_x),Ysz,Xsz);
		get_pol_center(cent_y, cent_x, fm->field[f].cg_y, fm->field[f].cg_x, &(fm->field[f].cg_ang), &(fm->field[f].cg_rad));
		fm->field[f].cg_rate = fm->rate_map[(int)fm->field[f].cg_y][(int)fm->field[f].cg_y];
		/* FIELD COHERENCE NOT IMPLEMENTED YET */
		fm->field[f].coh = -1.0;
	}
	
	if (fm->if_time)
		fm->if_rate = (float) fm->if_spks / fm->if_time;
	else 
		fm->if_rate = -1.0;
	if (fm->oof_time)
		fm->oof_rate = (float) fm->oof_spks / fm->oof_time;
	else 
		fm->oof_rate = -1.0;
	
	fm->in2out = fm->if_rate / fm->oof_rate;

	return;
}


get_field_center(fld, fnr, fmap, tm, sp,Ysz,Xsz,cent_y,cent_x)
FIELD	*fld;
int	fnr;
INT_MAP	fmap;
INT_MAP	tm;
INT_MAP	sp;
int	Ysz, Xsz;
float	cent_y, cent_x;
{
	int	i, j;
	float	yc_tmp, xc_tmp;
	int	ttot, stot;
	float	tmp, max_rt_tmp = 0.0;
	float	get_2x2_rate();
	void	get_pol_center();

	ttot = stot = 0;

	for (i = 0; i < Ysz -1; i++)
		for (j = 0; j < Xsz -1; j++)
			if (fmap[i][j] == fnr) {
				ttot += tm[i][j];
				stot += sp[i][j];
				if ((tmp = get_2x2_rate(i, j, Ysz, Xsz, fnr, fmap, tm, sp)) > max_rt_tmp) {
					max_rt_tmp = tmp;
					yc_tmp = (float)i +  0.5;
					xc_tmp = (float)j + 0.5;
				}
			}

	tmp = (float) stot / (float) ttot;
	fld->rate   = tmp * (float) PRAHA_SAMPS_PER_SEC;
	fld->cent_rate = max_rt_tmp * (float) PRAHA_SAMPS_PER_SEC;
	fld->cent_y = yc_tmp;
	fld->cent_x = xc_tmp;
/*
	fld->cent_y = cent_y - yc_tmp;
	fld->cent_x  = xc_tmp - cent_x;
*/

	get_pol_center(cent_y, cent_x, (float) yc_tmp, (float) xc_tmp, &(fld->cent_ang), &(fld->cent_rad));

}


void	get_pol_center(a_cent_y, a_cent_x, f_cent_y, f_cent_x, ang, rad)
float	a_cent_x, a_cent_y, f_cent_y, f_cent_x, *ang, *rad;
{
	double	hypot(), atan2(), dx, dy;
	
	if((f_cent_y) && (f_cent_x == 0.0)){
		*rad = 0.0;
		*ang = 0.0;
	}else{
		*rad = (float) hypot((double) (f_cent_x - a_cent_x), (double) (f_cent_y - a_cent_y));
		*ang = (float) atan2((double) (a_cent_y - f_cent_y), (double) (f_cent_x - a_cent_x));

		if (*ang < 0)
			*ang += TWO_PI;

		*ang *= DEG_PER_RAD;
	}
	return;
}
float   get_2x2_rate(y, x, Ysz, Xsz, f, fmap, tm, sp)
int     y, x, f;
INT_MAP	fmap;
INT_MAP	tm;
INT_MAP	sp;
int     Ysz, Xsz;
{
        int     tsum, ssum, i, j;
        float   tmp_rt = -1.0;

        tsum = ssum = 0;

        for (i = y;i < (y + 2); i++){
                for (j = x; j < (x + 2); j++){
                        if (f == fmap[i][j]) {
                                tsum += tm[i][j];
                                ssum += sp[i][j];
                        }
		}
	}
        tmp_rt = (float) ssum / (float) tsum;
 
        return tmp_rt;
}
