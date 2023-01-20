/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright (C) Martin Koehler, 2017-2023
 */

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <ctype.h>
#include "cscutils/hdf.h"

#include "benchmark.h"


int benchmark_sylv_setup_to_type(char *TRANSA, char*TRANSB, double sign)
{
    char ta = tolower(TRANSA[0]); 
    char tb = tolower(TRANSB[0]); 
    if ( sign == 1 ) {
        if ( ta == 'n' && tb == 'n' ) {
            return 0; 
        } else if ( ta == 'n' && tb == 't') {
            return 2; 
        } else if ( ta == 't' && tb == 'n') {
            return 4; 
        } else {
            return 6; 
        }
         
    } else {
        if ( ta == 'n' && tb == 'n' ) {
            return 1; 
        } else if ( ta == 'n' && tb == 't') {
            return 3; 
        } else if ( ta == 't' && tb == 'n') {
            return 5; 
        } else {
            return 7; 
        }
    }

}


int benchmark_gsylv_setup_to_type(char *TRANSA, char*TRANSB, double sign)
{
    char ta = tolower(TRANSA[0]); 
    char tb = tolower(TRANSB[0]); 
    if ( sign == 1 ) {
        if ( ta == 'n' && tb == 'n' ) {
            return 0; 
        } else if ( ta == 'n' && tb == 't') {
            return 2; 
        } else if ( ta == 't' && tb == 'n') {
            return 4; 
        } else {
            return 6; 
        }
         
    } else {
        if ( ta == 'n' && tb == 'n' ) {
            return 1; 
        } else if ( ta == 'n' && tb == 't') {
            return 3; 
        } else if ( ta == 't' && tb == 'n') {
            return 5; 
        } else {
            return 7; 
        }
    }

}

int benchmark_glyap_setup_to_type(char *TRANSA)
{
    char ta = tolower(TRANSA[0]); 
    if ( ta == 'n' ) 
        return 0; 
    else return 1; 

}

int benchmark_gstein_setup_to_type(char *TRANSA)
{
    char ta = tolower(TRANSA[0]); 
    if ( ta == 'n' ) 
        return 0; 
    else return 1; 

}


int benchmark_ggcsylv_setup_to_type(char *TRANSA, char*TRANSB, double sign)
{
    char ta = tolower(TRANSA[0]); 
    char tb = tolower(TRANSB[0]); 
    if ( sign == 1 ) {
        if ( ta == 'n' && tb == 'n' ) {
            return 0; 
        } else if ( ta == 'n' && tb == 't') {
            return 2; 
        } else if ( ta == 't' && tb == 'n') {
            return 4; 
        } else {
            return 6; 
        }
         
    } else {
        if ( ta == 'n' && tb == 'n' ) {
            return 1; 
        } else if ( ta == 'n' && tb == 't') {
            return 3; 
        } else if ( ta == 't' && tb == 'n') {
            return 5; 
        } else {
            return 7; 
        }
    }


}
