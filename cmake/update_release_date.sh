#!/bin/bash

DA=$(date +"%Y-%m-%d")
sed -i -e "s/^release-date:.*$/release-date: $DA/g" CODE
sed -i -e "s/^date-released: .*$/date-released: $DA/g" CITATION.cff

