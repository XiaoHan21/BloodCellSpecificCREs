#!/usr/bin/sh
infile=$1
outfile=$2

sed 's/\t/\n/g' $infile > $outfile && rm -rf ${infile}
