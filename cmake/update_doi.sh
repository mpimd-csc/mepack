#!/bin/bash
set -x

if [ $# -ne 1 ]; then
    echo "usage: $0 NEWDOI"
    exit 1
fi

NEW=$1

sed -i -e "s#^doi:.*\$#doi: $NEW#g" CITATION.cff
sed -i -e "s#^id:.*\$#id: $NEW#g" CODE
sed -i -e "s#DOI: .*\$#DOI: $NEW#g" README.md


