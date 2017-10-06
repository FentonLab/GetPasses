#include	<stdio.h>
#include	<stdlib.h>
#include	"tsf.h"

void	write_unit_sect(ofp, file_name, unit_num, fm)
	FILE	*ofp;
	FRAME_MAP	fm;
	char	*file_name;
	int	unit_num;
{
	char	*strcpy(), *strcat();
	
	fprintf(ofp,"\nSession: %s\t Unit %d\t Frame %d (Map number %d in tsf file)\n",file_name, unit_num, fm.map, fm.num);
	fprintf(ofp,"Arena: Center (%2.1f,%2.1f) Top %d Bot %d Left %d Right %d\n",fm.cent_y, fm.cent_x, fm.top_y, fm.bot_y, fm.lft_x, fm.rgt_x);
	fprintf(ofp,"%16s%10d%28s%10.2f\n","Area sampled:",fm.area,"Session Duration (s):",fm.tot_time);
	fprintf(ofp,"%16s%10d%28s%10.2f\n\n","Total Spikes:",fm.tot_spks,"Overall rate:",fm.grand_rate);
	fprintf(ofp,"%13s%7.3f%12s%7d\n","T-rate corr:",fm.rate_time_corr,"n = ",fm.n_rate_time_corr);
	fprintf(ofp,"%13s%7d%12s%7d%12s%7.4f%12s%7.4f%12s%7.3f\n","active pix:",fm.act_pix,"patches:",fm.patches,"coherence:",fm.z_coher,"info cnt:",fm.info_cnt, "z(ic):", fm.z_info);
	fprintf(ofp,"%13s%7.4f%12s%7.4f%12s%7.4f%12s%7.4f\n","active pix/A:",(float)fm.act_pix / (float)fm.area,"patches/A:",(float)fm.patches / (float)fm.area,"coherence/A:",fm.z_coher / (float)fm.area,"info rate:",fm.info_cnt * fm.grand_rate);
	
	return;
}

void	write_field_sect(ofp, fm)
	FILE	*ofp;
	FRAME_MAP	fm;
{
	char	*strcat();
	int	k;
        str256  temp, FldNum,FldRate,Size,SizeA,CentY,CentX,CentA,CentR;
        str256  FldTime,FldSpks,CGY,CGX,CGA,CGR,CGRate,FldCon,FldS2N,CentRate;

	fprintf(ofp,"Number of fields  = %d\n",fm.n_fields);
        	if(!fm.n_fields)
                	return;

	sprintf(FldNum,"\n%16s","Field #");
	sprintf(FldTime,"%16s","Field Time");
	sprintf(FldSpks,"%16s","Spikes in Field");
	sprintf(FldRate,"%16s","Field Rate");
	sprintf(FldS2N,"%16s","In rate:Out rate");
	sprintf(Size,"%16s","Size");
	sprintf(SizeA,"%16s","Size/Arena");
	sprintf(CentY,"%16s","Center Y");
	sprintf(CentX,"%16s","Center X");
	sprintf(CentA,"%16s","Center ang");
	sprintf(CentR,"%16s","Center rad");
	sprintf(CentRate,"%16s","Center Rate");
	sprintf(CGY,"%16s","Centroid Y");
	sprintf(CGX,"%16s","Centroid X");
	sprintf(CGA,"%16s","Centroid ang");
	sprintf(CGR,"%16s","Centroid rad");
	sprintf(CGRate,"%16s","Centroid Rate");
	sprintf(FldCon,"%16s","Concentration");
	sprintf(FldS2N,"%16s","Signal : Noise");

	fprintf(ofp,"In field:\tspikes %d\trate %2.1f\tin rate:out rate = %2.2f\n",fm.if_spks,fm.if_rate,fm.in2out);
	fprintf(ofp,"Out of field:\tspikes %d\trate %2.1f\n",fm.oof_spks,fm.oof_rate);

	for(k=0;k < fm.n_fields;k++){
		sprintf(temp,"%9d",k + 1);
		strcat(FldNum,temp);
		sprintf(temp,"%9.2f",fm.field[k].time);
		strcat(FldTime,temp);
		sprintf(FldSpks,"%9.2f",fm.field[k].spks);
		strcat(FldSpks,temp);
		sprintf(temp,"%9.2f",fm.field[k].rate);
		strcat(FldRate,temp);
		sprintf(temp,"%9.2f",(float)fm.field[k].s2n);
		strcat(FldS2N,temp);
		sprintf(temp,"%9d",fm.field[k].size);
		strcat(Size,temp);
		sprintf(temp,"%9.2f",(float)fm.field[k].size / (float)fm.area);
		strcat(SizeA,temp);
		sprintf(temp,"%9.2f",(float)fm.field[k].cent_y);
		strcat(CentY,temp);
		sprintf(temp,"%9.2f",(float)fm.field[k].cent_x);
		strcat(CentX,temp);
		sprintf(temp,"%9.2f",(float)fm.field[k].cent_ang);
		strcat(CentA,temp);
		sprintf(temp,"%9.2f",(float)fm.field[k].cent_rad);
		strcat(CentR,temp);
		sprintf(temp,"%9.2f",(float)fm.field[k].cent_rate);
		strcat(CentRate,temp);
		sprintf(temp,"%9.2f",(float)fm.field[k].cg_y);
		strcat(CGY,temp);
		sprintf(temp,"%9.2f",(float)fm.field[k].cg_x);
		strcat(CGX,temp);
		sprintf(temp,"%9.2f",(float)fm.field[k].cg_ang);
		strcat(CGA,temp);
		sprintf(temp,"%9.2f",(float)fm.field[k].cg_rad);
		strcat(CGR,temp);
		sprintf(temp,"%9.2f",(float)fm.field[k].cg_rate);
		strcat(CGRate,temp);
		sprintf(temp,"%9.2f",(float)fm.field[k].concent);
		strcat(FldCon,temp);
	}
	strcpy(temp,"\n");
	strcat(FldNum,temp);
	strcat(FldTime,temp);
	strcat(FldSpks,temp);
	strcat(FldRate,temp);
	strcat(FldS2N,temp);
	strcat(Size,temp);
	strcat(SizeA,temp);
	strcat(CentY,temp);
	strcat(CentX,temp);
	strcat(CentA,temp);
	strcat(CentR,temp);
	strcat(CentRate,temp);
	strcat(CGY,temp);
	strcat(CGX,temp);
	strcat(CGA,temp);
	strcat(CGR,temp);
	strcat(CGRate,temp);
	strcat(FldCon,temp);

	fprintf(ofp,"%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",FldNum,FldRate,FldS2N,Size,SizeA,CentY,CentX,CentA,CentR,CentRate,CGY,CGX,CGA,CGR,CGRate,FldCon,temp);
				
	return;
}

void	write_table(ofp,file_name,unit_num, fm)
	FILE	*ofp;
	char	*file_name;
	int	unit_num;
	FRAME_MAP	fm;
{
	int	k;

	fprintf(ofp,"%s\t%d\t%d\t",file_name, unit_num, fm.map);
		 /* session unit# frame# */
	fprintf(ofp,"%d\t%0.2f\t",fm.tot_spks,fm.grand_rate);
		/* #spikes mean_rate */
	fprintf(ofp,"%0.2f\t%d\t%0.2f\t",(float)fm.act_pix / (float)fm.area, fm.patches, fm.z_coher);
		/* p_active patches coherence */
	fprintf(ofp,"%0.2f\t%0.2f\t%0.2f\t", fm.info_cnt, fm.z_info, fm.info_cnt * fm.grand_rate); 
		/* info_cnt z(ic) info_rate */
	fprintf(ofp,"%0.2f\t%0.2f\t%0.2f\t%0.2f\t",fm.threshold, fm.if_rate, fm.oof_rate, fm.in2out);
		/* threshold infield_rate out_of_field_rate ratio */
	
	/* Field level analysis	*/
	for(k=0;k < fm.n_fields;k++){
		fprintf(ofp,"%s\t%d\t%d\t%d\t%0.2f\t",file_name, unit_num, fm.map, k + 1, fm.threshold);
			 /* session unit# map# field# threshold */
		fprintf(ofp,"%d\t%0.2f\t%0.2f\t",fm.field[k].spks, fm.field[k].rate, (float)fm.field[k].s2n);
			/* field spikes rate sig2noise */
		fprintf(ofp,"%0.2f\t",fm.field[k].concent);
			/* concentration */
		fprintf(ofp,"%d\t%0.2f\t",fm.field[k].size, (float)fm.field[k].size / (float)fm.area);
			/* size p_of_arena */
		fprintf(ofp,"%0.2f\t%0.2f\t",(float)fm.field[k].cg_y,(float)fm.field[k].cg_x);
			/* centroid y x */
		fprintf(ofp,"%0.2f\t%0.2f\t",(float)fm.field[k].cg_ang,(float)fm.field[k].cg_rad);
			/* centroid ang rad */
		fprintf(ofp,"%0.2f\t",(float)fm.field[k].cent_rate);
			/* center rate */
	}
	fprintf(ofp,"\n\n");
	return;
}
