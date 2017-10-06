#include <stdio.h>
#include "tsf.h"

count_one_color(map, Ysz, Xsz, color)
INT_MAP map;
int	Ysz, Xsz, color;

{
	void	make_color_map();
        INT_MAP     one_col_map;

        make_color_map(map,one_col_map,Ysz, Xsz, color);
        return(get_island_count(one_col_map, Ysz, Xsz));
}
void make_color_map(m1,m2,Ysz, Xsz, n)
INT_MAP     m1,m2;
int     Ysz, Xsz, n;

{
        int     i,j;

        for(i=0;i<Ysz;i++)
                for(j=0;j<Xsz;j++)
			m2[i][j] = 0;

        for(i=0;i<Ysz;i++){
                for(j=0;j<Xsz;j++){
                        if(m1[i][j] == n)
                                m2[i][j] = 1;
		}
	}
	return;
}

get_island_count(m, Ysz, Xsz)
INT_MAP     m;
int	Ysz, Xsz;
{
        int     ys,yf,xf,n_fields, first(), take_out_field();

        ys = n_fields = 0;

        while(1){
                first(Ysz, Xsz, ys,1,m,&yf,&xf);
                if(yf == -1) break;
                m[yf][xf] = 2;
                take_out_field(m,Ysz,Xsz);
                n_fields++;
                ys = yf;
                }
        return(n_fields);
}

/* static */ 
int	first(Ysz, Xsz, ys,seek,map,yf,xf)
int     Ysz, Xsz, ys,seek,*yf,*xf;
INT_MAP     map;

{
        register        i,j;

        *yf = -1;

        for(i=ys;i<Ysz;i++)
                for(j=0;j<Xsz;j++)
                        if(map[i][j] == seek){
                                *yf = i;
                                *xf = j;
                                return;
                                }
}

/* static */
int	 take_out_field(m,Ysz,Xsz)
INT_MAP     m;
int	Ysz,Xsz;

{
        int     y,x, first(), set_neighs();

        while(1){
                first(Ysz,Xsz,0,2,m,&y,&x);
                if(y == -1) break;
                m[y][x] = 3;
                set_neighs(m,y,x);
                }

}

/* static */
int	 set_neighs(map,y,x)
INT_MAP     map;
register        y,x;

{
        if(map[y-1][x] == 1) map[y-1][x] = 2;
        if(map[y][x-1] == 1) map[y][x-1] = 2;
        if(map[y][x+1] == 1) map[y][x+1] = 2;
        if(map[y+1][x] == 1) map[y+1][x] = 2;
}
