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

#include "source.h"

int source_subfix(char source_file[256]){
    /* check which kind of the source file format */
    char *pdot;
    int flag;
    
    pdot = strrchr(source_file, '.');
    if((pdot[1]=='s' && pdot[2]=='u')||(pdot[1]=='S' && pdot[2]=='U')){
        /* reads SU-format file */
        flag = 1;
    }else if(pdot[1]=='b' && pdot[2]=='i' && pdot[3]=='n'){
        /* reads BIN-format file */
        flag = 2;
    }else{
        flag = 0;
    }
    return flag;
}