#!/bin/bash

#
# Update the Copyright in all files.
#
set -x

X=$(pwd)

for d in "src" "examples" "matlab" "test"
do
cd "$d"

echo "Update Fortran"
grep -lRi "\![[:space:]]*copyright" | \
    xargs -n 1 \
    sed -i -e 's/![ \t]*Copyright.*$/! Copyright (C) Martin Koehler, 2017-2022/g'
echo "Update C"
grep -lRi " \* copyright" | \
    xargs -n 1 \
    sed -i -e 's/ \* Copyright.*$/ * Copyright (C) Martin Koehler, 2017-2022/g'
echo "Update Lua"
grep -lRi "\-- copyright"| \
    xargs -n 1  \
    sed -i -e 's/^-- Copyright.*$/-- Copyright (C) Martin Koehler, 2017-2022/g'
echo "Update MATLAB"
grep -lRi "% copyright"| \
    xargs -n 1  \
    sed -i -e 's/^% Copyright.*$/% Copyright (C) Martin Koehler, 2017-2022/g'


cd "$X"
done

sed -i -e 's/^Copyright [0-9-]* by/Copyright 2017-2022 by/g' README.md
