#!/bin/bash

if [ $# -ne 1 ]; then
    echo "usage: $0 NEWVERSION"
    exit 1
fi

NEW=$1

sed -i -e "s/^SET(VERSION.*)/SET(VERSION $NEW)/g" CMakeLists.txt
sed -i -e "s/^version:.*$/version: $NEW/g" CODE CITATION.cff
sed -i -e "s/^\*\*Version:\*\*.*$/**Version:** $NEW/g" README.md



