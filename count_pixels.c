#include "tsf.h"

int	count_gr0_pixels(map, Ysz, Xsz)
	INT_MAP	map;
	int	Ysz, Xsz;
{
	int	y, x, n =0;
	
	for(y=0; y < Ysz; y++)
		for(x=0; x < Xsz; x++)
			if(map[y][x] > 0)
				n++;
	return (n);
}
	
int	sum_gr0_pixels(map, Ysz, Xsz)
	INT_MAP	map;
	int	Ysz, Xsz;
{
	int	y, x, s =0;
	
	for(y=0; y < Ysz; y++)
		for(x=0; x < Xsz; x++)
			if(map[y][x] > 0)
				s += map[y][x];
	return (s);
}
void	zero_int(a,Y,X)
INT_MAP	a;
int	Y,X;
{
	int	y,x;

	for(y=0;y < Y;y++)
		for(x=0;x< X; x++)
			a[y][x] = 0;
	return;
}
void	copy_int(a,b,Y,X)
INT_MAP	a, b;
int	Y,X;
{
	int	y,x;

	for(y=0;y < Y;y++)
		for(x=0;x< X; x++)
			a[y][x] = b[y][x];
	return;
}
void	zero_float(a,Y,X)
FLOAT_MAP	a;
int	Y,X;
{
	int	y,x;

	for(y=0;y < Y;y++)
		for(x=0;x< X; x++)
			a[y][x] = 0.0;
	return;
}
void	copy_float(a,b,Y,X)
FLOAT_MAP	a, b;
int	Y,X;
{
	int	y,x;

	for(y=0;y < Y;y++)
		for(x=0;x< X; x++)
			a[y][x] = b[y][x];
	return;
}
