/*-------------------------------------------------------------------
 * Copyright (C) 2018  For the list of authors, see file AUTHORS.
 *
 * This file is part of SEISPLATFORM.
 * 
 * SEISPLATFORM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published 
 * by the Free Software Foundation, version 2.0 of the License only.
 * 
 * SEISPLATFORM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SEISPLATFORM. See file COPYING and/or 
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
-------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include "comm.h"

/* from https://blog.csdn.net/qq_41673920/article/details/81390972
 * function: replaces oldstr with newstr in [str] string
 * return: new [str]
 */
char *strrpc(char *str,char *oldstr,char *newstr){
    int i;
    char bstr[strlen(str)];//buffer
    memset(bstr,0,sizeof(bstr));
 
    for(i = 0; i < strlen(str); i++){
        if(!strncmp(str+i,oldstr,strlen(oldstr))){//search oldstr
            strcat(bstr,newstr);
            i += strlen(oldstr) - 1;
        }else{//save a char into buffer
            strncat(bstr,str + i,1);
        }
    }
 
    strcpy(str,bstr);
    return str;
}

void strrpc2(char str[256],char *oldstr,char *newstr, char str_out[256]){
    int i;
    char bstr[strlen(str)];//buffer
    memset(bstr,0,sizeof(bstr));
    memset(str_out,0,strlen(str_out));
 
    for(i = 0; i < strlen(str); i++){
        if(!strncmp(str+i,oldstr,strlen(oldstr))){//search oldstr
            strcat(bstr,newstr);
            i += strlen(oldstr) - 1;
        }else{//save a char into buffer
            strncat(bstr,str + i,1);
        }
    }
    strcpy(str_out, bstr);
}

/* copied from https://blog.csdn.net/xikangsoon/article/details/82054850 */
int StringToInteger(char *p)
{
	int value = 0;
	while (*p != '\0')
	{
		if ((*p >= '0') && (*p <= '9'))
		{
			value = value * 10 + *p - '0';
		}
		p++;
 
	}
	return value;
}

char **read_string(int *nrow, char filename[256]){
    
    int c;
    int nr, irow;
    char strtmp[256], **str_list = NULL;
    FILE *fp;
    
    fp=fopen(filename,"r");
    if (fp==NULL){
        printf(" File %s could not be opened !",filename);
        exit(0);
    }
    nr=0;
    /* The last line should be empty!!! */
    while ((c=fgetc(fp)) != EOF)
        if (c=='\n') ++(nr);
    rewind(fp);
    str_list = chmatrix(0, nr-1, 0, 255);
    for(irow=0;irow<nr;irow++){
        fscanf(fp,"%s %s",strtmp,str_list[irow]);
    }
    fclose(fp);
    printf(" Number of data found: %i\n",nr);
    *nrow= nr;
    
    return str_list;
}

char **read_string1(int *nrow, char filename[256]){
    
    int c;
    int nr, irow;
    char **str_list = NULL;
    FILE *fp;
    
    fp=fopen(filename,"r");
    if (fp==NULL){
        printf(" File %s could not be opened !",filename);
        exit(0);
    }
    nr=0;
    /* The last line should be empty!!! */
    while ((c=fgetc(fp)) != EOF)
        if (c=='\n') ++(nr);
    rewind(fp);
    str_list = chmatrix(0, nr-1, 0, 255);
    for(irow=0;irow<nr;irow++){
        fscanf(fp,"%s",str_list[irow]);
    }
    fclose(fp);
    printf(" Number of data found: %i\n",nr);
    *nrow= nr;
    
    return str_list;
}

/*
 gcc -o string_opts string_opts.c
 ./string_opts
 */
// int main(void){
//     char str[] = "Hello,China!";
//     char *str2=NULL;
//     char str_out[256];
//     printf("len_str=%d\n",(int) strlen(str));
//     strrpc(str,"China","World123");
//     printf("%s\n",str);
//     strrpc(str,"Hello","Hi");
//     printf("%s\n",str);
//     strrpc(str,"Hi,World","Hello,world");
//     printf("%s\n",str);
//     str2=strrpc(str,"Hello","HHi");
//     printf("%s\n",str2);
//     strrpc2(str2,"HHi","HHHHHH",str_out);
//     printf("\n%s\n",str2);
//     printf("%s\n",str_out);
//     
//     strrpc2("test_vx_proc123.bin","vx","psp",str_out);
//     printf("%s\n",str_out);
//     return 0;
// }

