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

#include "comm.h"

void update_shot_num(char filein[256],int ishot,char fileout[256]){
    /*updates the shot number in the path of source file*/
    char str_pre[256];
    char *pshot,*pdot;
    int num_pre;
    
    memset(str_pre, '\0', sizeof(str_pre)); 
    pshot = strstr(filein, "shot");
    if(pshot!=NULL){
        //pdot = strstr(filein, ".");
        pdot = strrchr(filein, '.');
        //printf("The substring is: %s\n", pshot);
        num_pre=pshot+4-filein;
        //printf("%d\n",num_pre);
        strncpy(str_pre,filein, num_pre);
        //printf("The substring is: %s\n",str_pre);
        sprintf(fileout,"%s%d%s",str_pre,ishot,pdot);
        //printf("The substring is: %s\n",fileout);
    }else{
        sprintf(fileout,"%s",filein);
    }
}

void update_shot_num2(char filein[256],int ishot,char fileout[256]){
    /*updates the shot number in the path of source file*/
    char str_pre[256];
    char *pshot,*pdot;
    int num_pre;
    
    memset(str_pre, '\0', sizeof(str_pre)); 
    pshot = strstr(filein, "shot");
    if(pshot!=NULL){
        //pdot = strstr(filein, ".");
        pdot = strchr(pshot, '.');
        //printf("The substring is: %s\n", pshot);
        num_pre=pshot+4-filein;
        //printf("%d\n",num_pre);
        strncpy(str_pre,filein, num_pre);
        //printf("The substring is: %s\n",str_pre);
        sprintf(fileout,"%s%d%s",str_pre,ishot,pdot);
        //printf("The substring is: %s\n",fileout);
    }else{
        sprintf(fileout,"%s",filein);
    }
}

int extract_shot_num(char filein[256]){
    /*updates the shot number in the path of source file*/
    char str_s[20];
    char *pshot,*pdot;
    int num_l, ishot;
    
    memset(str_s, '\0', sizeof(str_s)); 
    pshot = strstr(filein, "shot");
    if(pshot!=NULL){
        //pdot = strstr(filein, ".");
        pdot = strchr(pshot, '.');
        //printf("%s\n %s\n",pshot,pdot  );
        //printf("The substring is: %s\n", pshot);
        num_l=pdot-4-pshot;
        //printf("%d\n",num_pre);
        strncpy(str_s,pshot+4, num_l);
        //printf("%s\n",str_s);
        ishot = atoi(str_s);
    }else{
        printf("Can't extract the number from file %s",filein);
        exit(1);
    }
    return ishot;
}

// int main () {
//     //char SOURCE_FILE[256] = "/home/tao/IFOS_OPTIM/IFOS2D/src/toy_example_p_shot32.su";
//     char SOURCE_FILE[256] = "./snap/inc/./././test.vx.shot199.2.1.bin";
//     char str_pre[256];
//     char srcname[60],srcname2[60];
//     char *pshot,*pdot;
//     int num_pre,ishot=999;
// 
// //     memset(str_pre, '\0', sizeof(str_pre)); 
// //     pshot = strstr(SOURCE_FILE, "shot");
// //     if(pshot!=NULL){
// //         //pdot = strstr(SOURCE_FILE, ".");
// //         pdot = strrchr(SOURCE_FILE, '.');
// //         //printf("The substring is: %s\n", pshot);
// //         num_pre=pshot+4-SOURCE_FILE;
// //         printf("%d\n",num_pre);
// //         strncpy(str_pre,SOURCE_FILE, num_pre);
// //         printf("The substring is: %s\n",str_pre);
// //         sprintf(srcname,"%s%d%s",str_pre,ishot,pdot);
// //     }else{
// //         sprintf(srcname,"%s",SOURCE_FILE);
// //     }
// //     printf("The substring is: %s\n",srcname);
// //     
// //     update_shot_num(SOURCE_FILE,ishot,srcname2);
// //     printf("%s\n",srcname2);
//     
//     printf("%d\n",extract_shot_num(SOURCE_FILE));
//     return(0);
// }

