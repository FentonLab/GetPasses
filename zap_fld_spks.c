#include "map_types.h"

void	zap_out_of_fld_spks(spks, flds, fld_nr)
SPKS	spks;
PRT_MAP	flds;
int	fld_nr;
{
	int	i, j;

	for (i = 0; i < YSZ; i++)
		for (j = 0; j < XSZ; j++)
			if (flds[i][j] != fld_nr)
				spks[i][j] = 0;
	return;
}


void	zap_in_fld_spks(spks, flds, fld_nr)
SPKS	spks;
PRT_MAP	flds;
int	fld_nr;
{
	int	i, j;

	for (i = 0; i < YSZ; i++)
		for (j = 0; j < XSZ; j++)
			if (flds[i][j] == fld_nr)
				spks[i][j] = 0;
	return;
}


