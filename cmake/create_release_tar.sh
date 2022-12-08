#!/bin/bash
# VERSION=$(grep "^version:" CITATION.cff | cut -d' ' -f 2)
VERSION=$(git describe | sed 's/^v//g')
OUTPUT=../mepack-${VERSION}.tar.gz
git ls-files --recurse-submodules | \
    grep -v '.gitmodules' | \
    grep -v 'compile_commands' | \
    grep -v 'gitlab-ci' | \
    tar cazf "${OUTPUT}" --transform="s#^#mepack-${VERSION}/#g" -T-
echo "Wrote ${OUTPUT}"
