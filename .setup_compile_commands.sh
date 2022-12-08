#!/usr/bin/env bash
X=`mktemp -d cmake_tmp_XXXXXXXXXXXXXX`
echo "Generating compile commands in ${X}"
cmake -H. -B"${X}" -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=YES -DMATLAB=ON
cp "${X}/compile_commands.json" .
rm -rf "${X}"
