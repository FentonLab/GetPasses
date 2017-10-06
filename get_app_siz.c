#include	"bwmap.h"

get_app_size(m, n_regs, app_token)
APP_MAP	m;
int	n_regs, *app_token;

{
	int	i, j;
	int	region_sizes[AREA];
	int	max_size = 0;

	if (n_regs == 1) {
		for (i = 0; i < YSZ; i++)
			for (j = 0; j < XSZ; j++)
				if (m[i][j])
					max_size++;
		*app_token = 1;
		return(max_size);
	}

	for (i = 0; i < AREA; i++)
		region_sizes[i] = 0;

	for (i = 0; i < YSZ; i++)
		for (j = 0; j < XSZ; j++)
			if (m[i][j])
				region_sizes[m[i][j]]++;

	for (i = 1; i <= n_regs; i++)
		if (region_sizes[i] > max_size) {
			max_size = region_sizes[i];
			*app_token = i;
		}

	return max_size;
}


