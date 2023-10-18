#!/usr/bin/env sh
aspell -c --add-filter-path=$(pwd)/cmake --mode=ccpp -l en -p $(pwd)/cmake/aspell.pws $1


