#!/bin/bash

find src/ -type f | xargs -n 1 sed -i -e "s/[@\]date.*$/\\\\date `date +'%B %Y'`/g"
