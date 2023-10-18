#!/bin/bash
aspell -c --add-filter-path=$(pwd)/cmake --mode=fortran -l en -p $(pwd)/cmake/aspell.pws $1
