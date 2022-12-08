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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "FCMangle.h"
#include "mepack.h"
#include "mepack_internal.h"
#include "cscutils/lua.h"


static int mepack_lua_init_state = 0;
static csc_lua_t lua = NULL;
/**
 * \internal
 * \brief finalize the lua interpreter
 */
static void mepack_lua_exit(void)
{
    if (mepack_lua_init_state == 0)
        return;
    csc_lua_finalize(lua);
    mepack_lua_init_state = 0;
    return;
}


extern const char * mepack_default_config(void);

/**
 * \internal
 * \brief Init the LUA interface from Fortran
 */
void FC_GLOBAL_(mepack_lua_init,MEPACK_LUA_INIT)(Int *_info)
{
    char *luaconfig;
    FILE *fp;

    *_info = 0;

    if ( mepack_lua_init_state ) return;

    luaconfig = getenv("MEPACK_LUA_CONFIG");
    if ( luaconfig ) {

        fp = fopen(luaconfig, "r");
        if (!fp ) {
            fprintf(stderr, "MEPACK_LUA_CONFIG: %s is not readable. Settings are ignored.\n", luaconfig);
            mepack_lua_init_state = 0;
            *_info = 1;
            return;
        }
        fclose(fp);

        lua = csc_lua_init();
        mepack_lua_init_state = 1;

        if ( csc_lua_loadfile(lua, luaconfig) ) {
            fprintf(stderr, "MEPACK_LUA_CONFIG: Failed to load %s. Check its syntax.\n", luaconfig);
            *_info = 2;
            mepack_lua_exit();
            return;
        }
        if ( mepack_verbose_get()) {
            fprintf(stderr, "MEPACK: Loaded configuration from %s\n", luaconfig);
        }
    } else {
        if ( mepack_verbose_get()) {
            /* Load the default configuration */
            fprintf(stderr, "MEPACK: MEPACK_LUA_CONFIG not set. Using defaults.\n");
        }
        lua = csc_lua_init();
        mepack_lua_init_state = 1;

        if ( csc_lua_loadstring(lua, mepack_default_config()) ) {
            fprintf(stderr, "MEPACK: Loading default configuration failed.\n");
            *_info = 2;
            mepack_lua_exit();
            return;
        }
    }

    if ( csc_lua_run(lua)) {
        fprintf(stderr, "MEPACK_LUA_CONFIG: Failed to execute %s. Check its syntax.\n", luaconfig);
        *_info = 3;
        mepack_lua_exit();
        return;
    }

    mepack_lua_init_state = 1;
    atexit(mepack_lua_exit);
    return;
}

/**
 * \internal
 *
 * Return a pointer to the state of the LUA interpreter.
 */
void *mepack_lua_state(void) {
    return (void*) lua;
}

/**
 * \internal
 *
 * Return a Boolean value if the LUA interpreter is loaded.
 */
int mepack_lua_loaded(void)
{
    return (int) mepack_lua_init_state;
}
