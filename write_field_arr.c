#include <stdio.h>
#include <stdlib.h>
#include "tsf.h"


void	write_field_Iarr(file, arr, field_map, field, Ysz, Xsz)
char	*file;
INT_MAP	arr, field_map;
int	field, Ysz, Xsz;
{
	FILE	*fp, *fopen();
	int	i, j;

	if((fp = fopen(file, "w")) == NULL){       
        	fprintf(stderr,"Can't open array output file: %s\n",file);
               	exit(1);
        }

	for(i = 0; i < Ysz; i++){
		for(j = 0; j < Xsz; j++){
			if(field_map[i][j] <= 0){
				fprintf(fp, "%d\n", field_map[i][j]);
				continue;
			}
			if(field_map[i][j] == field){
				fprintf(fp, "%d\n", arr[i][j]);
			}else{
				fprintf(fp, "%d\n", 0);
			}
		}
	}
	
	fclose(fp);
	return;
}

void	write_field_Farr(file, arr, field_map, field, Ysz, Xsz)
char	*file;
FLOAT_MAP	arr;
INT_MAP	field_map;
int	field, Ysz, Xsz;
{
	FILE	*fp, *fopen();
	int	i, j;

	if((fp = fopen(file, "w")) == NULL){       
        	fprintf(stderr,"Can't open array output file: %s\n",file);
               	exit(1);
        }

	for(i = 0; i < Ysz; i++){
		for(j = 0; j < Xsz; j++){
			if(field_map[i][j] <= 0){
				fprintf(fp, "%0.2f\n",(float)field_map[i][j]);
				continue;
			}
			if(field_map[i][j] == field){
				fprintf(fp, "%0.2f\n",arr[i][j]);
			}else{
				fprintf(fp, "%0.2f\n", 0.0);
			}
		}
	}
	
	fclose(fp);
	return;
}

void	write_Iarr(file, arr, Ysz, Xsz)
char	*file;
INT_MAP	arr;
int	Ysz, Xsz;
{
	FILE	*fp, *fopen();
	int	i, j;

	if((fp = fopen(file, "w")) == NULL){       
        	fprintf(stderr,"Can't open array output file: %s\n",file);
               	exit(1);
        }

	for(i = 0; i < Ysz; i++)
		for(j = 0; j < Xsz; j++)
			fprintf(fp, "%d\n", arr[i][j]);
	
	fclose(fp);
	return;
}
