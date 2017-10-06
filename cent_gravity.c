#include        "tsf.h"

void	cent_gravity(r,f,field_nr,y,x, Ysz,Xsz)
FLOAT_MAP	r;
INT_MAP		f;
int	field_nr, Ysz,Xsz;
float	*y,*x;
{
        int     i, j, n = 0;
	float	mean_rate = 0.0;
	
	*x = *y = 0.0;

        for (i = 0; i < Ysz; i++)
                for (j = 0; j < Xsz; j++)
                        if (f[i][j] == field_nr){
				*y += (float)i * r[i][j];
				*x += (float)j * r[i][j];
				mean_rate += r[i][j];
				n++;
                        }

	*y /= (float) n;
	*x /= (float) n;
	mean_rate /= (float)n;

	*y /= mean_rate;
	*x /= mean_rate;

	*y += 0.5;
	*x += 0.5;
	return;
}
