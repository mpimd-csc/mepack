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
 * Copyright (C) Martin Koehler, 2017-2022
 */
/* For Documentation see xerror_handler_doc.f90  */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h> 
#include <string.h>

#ifndef FC_GLOBAL
#include "FCMangle.h"
#endif 

#ifndef Int 
#define Int int 
#endif 

static void (*error_handler_c)(const char *, int ) = NULL; 
static void (*error_handler_f)(const char *, Int* , Int ) = NULL; 
void FC_GLOBAL_(xerror_handler,XERROR_HANDLER)(const char * str, Int * info, Int length) {
    char *istr; 

    istr = malloc(sizeof(char)*(length+1)); 
    memcpy(istr, str, length*sizeof(char)); 
    istr[length] = '\0'; 
    if ( error_handler_c != NULL ) {
        int _info = (int) *info; 
        error_handler_c(istr, _info); 
    } else if ( error_handler_f != NULL ) {
        error_handler_f(str, info, length); 
    } else {
        fprintf(stderr, "On entry function %s has an invalid argument %d.\n", istr, *info);
    }


    free(istr); 
    return; 
}

void FC_GLOBAL_(xerror_set_handler_c,XERROR_SET_HANDLER_C) ( void (*eh)(const char*,int)) 
{
    error_handler_f = NULL; 
    error_handler_c = eh; 
    return; 
}

void FC_GLOBAL_(xerror_set_handler_f,XERROR_SET_HANDLER_F) ( void (*eh)(const char*,Int*, Int)) 
{
    error_handler_f = eh; 
    error_handler_c = NULL; 
    return; 
}





