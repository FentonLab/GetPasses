/*	routines for reading Kaminsji style data files	*/

#include <stdio.h>
#include "tsf.h"

char    *get_string(fp, d1, d2) /* gets text starting and ending with d1 and d2 */
        FILE    *fp;
        BYTE    d1, d2;
{
        static char     s[256];
        int     i;
        BYTE    c;

        i = 0;
        while((c = (BYTE)getc(fp)) != EOF)
                if(c == d1){
                	sprintf(s,"%c",c);
			i++;
                        break;
		}
        while((c = (BYTE)getc(fp)) != EOF){
                sprintf(s+i,"%c",c);
                if(c == d2){
			i++;
                        break;
		}
                i++;
	}
        sprintf(s+i,"%c",'\0');
        return s;
}

char    *get_delimited_string(fp, d1, d2) /* gets text between d1 and d2 */
        FILE    *fp;
        BYTE    d1, d2;
{
        static char     s[256];
        int     i;
        BYTE    c;

        while((c = (BYTE)getc(fp)) != EOF)
                if(c == d1)
                        break;
        i = 0;
        while((c = (BYTE)getc(fp)) != EOF){
                if(c == d2)
                        break;
                sprintf(s+i,"%c",c);
                i++;
        }
        sprintf(s+i,"%c",'\0');
        return s;
}

BYTE	check_comma(fp)
	FILE *fp;
{	

	return((BYTE)(getc(fp) == ','));
}

int	ssub_char(s,old,new)
	char	*s,old, new;
{
	int	i = 0, n= 0;

	while(*(s+i)!= '\0'){
		if(*(s+i) == old){
			*(s+i) = new;
			n++;
		}
		i++;
	}
	return(n);
}
int	fsub_char(fp,old,new)
	FILE	*fp;
	char	old, new;
{
	char	c;
	int	i=0;

	while((c = (char)getc(fp)) != EOF){
		if(c == old){
			fseek(fp,-1L,1);	/* go back one */
			putc(new, fp);
			fseek(fp,1L,1);		/* go forward one */
			i++;
		}
	}
	return (i);
}

char	*trunc_str(s,c)
	char *s;
	int	c;
{
	while(*s != c)
		s++;
	*s = '\0';
	return s;
}
char	*get_line(f,s)
	FILE	*f;
	char *s;
{
	int	i=0;

	while(fscanf(f,"%c",s+i) != EOF){
		if(*(s+i) == '\n'){
			*(s+i+1) = '\0';
			return s;
		}
		i++;
	}
	return NULL;
}
