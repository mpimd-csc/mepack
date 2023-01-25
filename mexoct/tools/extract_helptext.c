/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
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
 * Copyright (C) Martin Koehler, 2019
 */



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <dlfcn.h>
#include "mex.h"

#  ifdef __cplusplus
#    define MEXOCT_EXTERN_C extern "C"
#  else
#    define MEXOCT_EXTERN_C extern
#  endif
#ifdef _MSC_VER
#define MEXOCT_DLL_EXPORT_SYM __declspec(dllexport)
#define MEXOCT_DLL_IMPORT_SYM __declspec(dllimport)
#elif __GNUC__ >= 4
#define MEXOCT_DLL_EXPORT_SYM __attribute__((visibility("default")))
#define MEXOCT_DLL_IMPORT_SYM __attribute__((visibility("default")))
#else
#define MEXOCT_DLL_EXPORT_SYM
#define MEXOCT_DLL_IMPORT_SYM
#endif

MEXOCT_EXTERN_C MEXOCT_DLL_EXPORT_SYM
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2 ) {
        mexErrMsgTxt("two arguments required.");
    }

    char *infi = mxArrayToString(prhs[0]);
    char *outfi = mxArrayToString(prhs[1]);
    dlerror();
    void *handle = dlopen ( infi, RTLD_LAZY);
    if ( !handle) {
        mexPrintf("Failed to open %s ( error = %s ) \n", infi, dlerror());
        return;
    }


    dlerror();
    typedef char * (*gethelpstring_t)(void);
    gethelpstring_t gethelpstring;
    *(void **) &gethelpstring =  dlsym(handle, "mexoct_gethelptext");
    if ( !gethelpstring){
        mexPrintf("No help string found in %s ( error = %s )\n", infi, dlerror());
        return;
    }
    char *helpstring = strdup(gethelpstring());
    /* mexPrintf("%s\n", helpstring); */

    FILE *out;
    out = fopen ( outfi, "w");
    if ( !out) {
        free(helpstring);
        mexPrintf("Failed to open %s for writing.\n", outfi);
        exit(4);
    }

    size_t len = strlen(helpstring);
    size_t i;
    fprintf(out, "%%%% ");
    for (i = 0; i < len; i++) {
        if ( helpstring[i] == '\n' ) {
            fprintf(out, "\n%%");
        } else {
            fprintf(out, "%c", helpstring[i]);
        }
    }

    free(helpstring);
    fclose(out);

    dlclose(handle);
    return ;

}

