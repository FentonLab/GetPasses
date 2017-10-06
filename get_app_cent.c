#include	"tsf.h"

find_app_center(fm, time, Ysz, Xsz)
FRAME_MAP	*fm;
INT_MAP	time;
int	Ysz,Xsz;
{
	int	i, j;

	for (i = 0; i < Ysz; i++)
		for (j = 0; j < Xsz; j++)
			if (time[i][j] > -1) 
				goto found_top;

found_top:
	;

	fm->top_y = i;

	for (i = Ysz - 1; i >= 0; i--)
		for (j = 0; j < Xsz; j++)
			if (time[i][j] > -1) 
				goto found_bot;

found_bot:
	;

	fm->bot_y = i;

	for (j = 0; j < Xsz; j++)
		for (i = 0; i < Ysz; i++)
			if (time[i][j] > -1) 
				goto found_lft;

found_lft:
	;

	fm->lft_x = j;

	for (j = Xsz - 1; j >= 0; j--)
		for (i = 0; i < Ysz; i++)
			if (time[i][j] > -1) 
				goto found_rgt;

found_rgt:
	;

	fm->rgt_x = j;

	fm->cent_y = ((double) fm->top_y + (double) fm->bot_y + 1.0) / 2.0;
	fm->cent_x = ((double) fm->lft_x + (double) fm->rgt_x + 1.0) / 2.0;
}


