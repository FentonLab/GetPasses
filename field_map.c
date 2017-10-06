#include  "bwmap.h"

static neighbors_1();
static int	compare_1();

get_field_map(rate, intmap, thr)
RATE	rate;
FLD_MAP	intmap;
float	thr;

{
	int	i, j, intval = 1;

	zero_usi(intmap);

	for (i = 1; i < YSZ - 1; i++)
		for (j = 1; j < XSZ - 1; j++)
			if (rate[i][j] >= thr)
				if (intmap[i][j] < 1) {
					neighbors_1(i, j, rate, intmap, thr, intval);
					intval++;
				}
	return(intval - 1);
}


neighbors_1(i, j, r, intmap, thr, intval)
int	i, j, intval;
RATE	r;
float	thr;
u_map	intmap;

{
	int	itmp, jtmp;

	intmap[i][j] = intval;
	itmp = i;

	for (i -= 1; i < itmp + 2; i += 2)
		if (intmap[i][j] < 1)
			if (r[i][j] >= thr)
				neighbors_1(i, j, r, intmap, thr, intval);

	i = itmp;
	jtmp = j;

	for (j -= 1; j < jtmp + 2; j += 2)
		if (intmap[i][j] < 1)
			if (r[i][j] >= thr)
				neighbors_1(i, j, r, intmap, thr, intval);
}


typedef struct {
	int	f_nr;
	int	f_sz;
} sort_flds;

sort_flds	*to_sort;

prune_field_list(fmap, n_flds, min_siz)
FLD_MAP	fmap;
int	n_flds, min_siz;
{
	void	 * calloc();
	int	i, j, n;
	int	tmp, good_flds;
	FLD_MAP	t_map;
	void	qsort();

	to_sort = (sort_flds * ) calloc((unsigned)n_flds, sizeof(sort_flds));

	for (n = 1; n <= n_flds; n++)
		to_sort[n-1].f_nr = n;

	for (i = 0; i < YSZ; i++)
		for (j = 0; j < XSZ; j++)
			if (fmap[i][j])
				to_sort[fmap[i][j] - 1].f_sz++;

	qsort((void * ) to_sort, (unsigned)n_flds, sizeof(sort_flds), compare_1);

	n = 0;
	while (n < n_flds && to_sort[n].f_sz >= min_siz) 
		n++;

	good_flds = n;

	/* This next loop need not be done if n == n_flds, but what the heck */

	for (i = 0; i < YSZ; i++)
		for (j = 0; j < XSZ; j++)
			if (tmp = fmap[i][j])
				for (n = good_flds; n < n_flds; n++)
					if (tmp == to_sort[n].f_nr) {
						fmap[i][j] = 0;
						break;
					}

	zero_usi(t_map);

	for (i = 0; i < YSZ; i++)
		for (j = 0; j < XSZ; j++) {
			for (n = 0; n < good_flds; n++)
				if (fmap[i][j] == to_sort[n].f_nr) {
					t_map[i][j] = n + 1;
					break;
				}
		}

	copy_usi(fmap, t_map);

	return(good_flds);
}


compare_1(a, b)
sort_flds	*a, *b;
{
	if (a->f_sz > b->f_sz) 
		return(-1);
	if (a->f_sz < b->f_sz) 
		return(1);
	return(0);
}


void	do_field_stuff(sess, c)
SESSION	sess;
CLUSTER	*c;
{
	void	 * calloc();
	int	f, fnr, i, j, oof_spks, oof_time;

	make_show_fld_map(sess, c);

	c->fields = (FIELD * ) calloc((unsigned) c->n_fields, sizeof(FIELD));

	oof_spks = 0;
	oof_time = 0;
	for (i = 0; i < YSZ; i++)
		for (j = 0; j < XSZ; j++)
			for (f = 0, fnr = 1; f < c->n_fields; f++, fnr++) {
				if (fnr == c->field_map[i][j]) {
					c->fields[f].size++;
					c->fields[f].f_tm += sess.time[i][j];
					c->fields[f].f_sp += c->spks[i][j];
				}
				else {
					oof_spks += c->spks[i][j];
					oof_time += sess.time[i][j];
				}
			}

	for (f = 0, fnr = 1; f < c->n_fields; f++, fnr++) {
		c->fields[f].f_rt = (float) c->fields[f].f_sp / (float)c->fields[f].f_tm;
		c->fields[f].f_rt *= (float) SAMP_PER_SEC;
		c->fields[f].field_nr = fnr;
		get_field_center(sess, &(c->fields[f]), fnr, c->field_map, sess.time, c->spks);
	}
	
	if (oof_time)
		c->out_of_field_rate = ((float) oof_spks * (float) SAMP_PER_SEC) / (float) oof_time;
	else 
		c->out_of_field_rate = 0.0;

	return;
}


get_field_center(s, fld, fnr, fmap, tm, sp)
SESSION	s;
FIELD	*fld;
int	fnr;
FLD_MAP	fmap;
TIME	tm;
SPKS	sp;
{
	int	i, j, yc_tmp, xc_tmp;
	int	ttot, stot;
	float	tmp, max_rt_tmp = 0.0;
	float	get_3x3_rate();
	float	ang, rad;

	ttot = stot = 0;

	for (i = 0; i < YSZ; i++)
		for (j = 0; j < XSZ; j++)
			if (fmap[i][j] == fnr) {
				ttot += tm[i][j];
				stot += sp[i][j];
				if ((tmp = get_3x3_rate(i, j, Ysz, Xsz, fnr, fmap, tm, sp)) > max_rt_tmp) {
					max_rt_tmp = tmp;
					yc_tmp = i;
					xc_tmp = j;
				}
			}

	tmp = (float) stot / (float) ttot;
	fld->f_rt   = tmp * (float) SAMP_PER_SEC;
	fld->c_rt   = max_rt_tmp * (float) SAMP_PER_SEC;
	fld->y_ctr  = yc_tmp;
	fld->x_ctr  = xc_tmp;

	get_pol_center(s.app_cent_y, s.app_cent_x, (float) yc_tmp, (float) xc_tmp, &ang, &rad);

	fld->ang_ctr = ang;
	fld->rad_ctr = rad;
}


float	get_3x3_rate(y, x, Ysz, Xsz, f, fmap, tm, sp)
int	y, x, f;
FLD_MAP	fmap;
TIME	tm;
SPKS	sp;
int	Ysz, Xsz;
{
	int	tsum, ssum, i, j;
	float	tmp_rt;

	tsum = ssum = 0;

	for (i = y - 1; i < y + 2; i++)
		for (j = x - 1; j < x + 2; j++)
			if (f == fmap[i][j]) {
				tsum += tm[i][j];
				ssum += sp[i][j];
			}

	tmp_rt = (float) ssum / (float) tsum;

	return tmp_rt;
}


get_pol_center(a_cent_y, a_cent_x, f_cent_y, f_cent_x, ang, rad)
float	a_cent_x, a_cent_y, f_cent_y, f_cent_x, *ang, *rad;
{
	double	hypot(), atan2();

	*rad = (float) hypot((double) (f_cent_x - a_cent_x), (double) (f_cent_y - a_cent_y));
	/*  fprintf(stderr,"%f         %f         %f          %f\n",a_cent_y,a_cent_x,f_cent_y,f_cent_x);  */
	*ang = (float) atan2((double) (a_cent_y - f_cent_y), (double) (f_cent_x - a_cent_x));

	if (*ang < 0)
		*ang += TWO_PI;

	*ang *= DEG_PER_RAD;
}


make_show_fld_map(sess, c)
SESSION	sess;
CLUSTER	*c;
{
	int	i, j;

	for (i = 0; i < YSZ; i++)
		for (j = 0; j < XSZ; j++)
			if (sess.app_image[i][j])
				c->show_flds[i][j] = c->field_map[i][j] + 1;
			else
				c->show_flds[i][j] = 0;
}


