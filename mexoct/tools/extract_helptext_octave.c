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
#include <errno.h>

#include <dlfcn.h>

int main(int argc, char **argv)
{
    char *infile = NULL;
    char *outfile = NULL;
    if ( argc != 3) {
        fprintf(stderr, "usage: %s input.oct output.m\n", argv[0]);
        return 1;
    }


    dlerror();
    infile = strdup(argv[1]);
    outfile = strdup(argv[2]);
    void *handle = dlopen ( infile, RTLD_LAZY | RTLD_LOCAL);
    if (!handle) {
        fprintf(stderr, "Failed to open %s - error: %s\n", infile, dlerror());
        return 1;
    }

    typedef char * (*gethelpstring_t)(void);
    gethelpstring_t help = NULL;

    *(void **) &help = dlsym(handle, "octHelp");
    if (!help) {
        fprintf(stderr, "Failed to get help text. error: %s\n", dlerror());
        dlclose(handle);
        return 1;
    }

    FILE *fpout = NULL;
    fpout = fopen(outfile, "w");
    if (!fpout) {
        int err = errno;
        fprintf(stderr, "Failed to open file %s. error: %s\n", outfile, strerror(err));
        dlclose(handle);
        return 1;
    }
    char *helpstring = strdup(help());
    size_t len = strlen(helpstring);
    size_t i;
    fprintf(fpout, "%%%% ");
    for (i = 0; i < len; i++) {
        if ( helpstring[i] == '\n' ) {
            fprintf(fpout, "\n%% ");
        } else {
            fprintf(fpout, "%c", helpstring[i]);
        }
    }

    free(helpstring);
    fclose(fpout);
    fprintf(stdout, "Help text from %s written to %s\n", infile, outfile);
    dlclose(handle);
    return 0;


}
