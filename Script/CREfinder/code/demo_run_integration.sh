#!/usr/bin/sh
## This is a demo for Software SpecificCREsIdentification
## Written by XiaoHan

echo "start at `date`" &&

/jdfsbjcas1/ST_BJ/P21H28400N0232/xiaohan2/Software/miniconda/envs/R4.2/bin/Rscript CREfinder.r \
--Mode Integration \
--PeakGene /jdfsbjcas1/ST_BJ/P21H28400N0232/xiaohan2/Software/CREfinder/Demo/01-input/Integration/i1.peak2gid2symbol.txt \
--AUCtbl /jdfsbjcas1/ST_BJ/P21H28400N0232/xiaohan2/Software/CREfinder/Demo/01-input/Integration/i2.celltype2auc.txt \
--Outdir /jdfsbjcas1/ST_BJ/P21H28400N0232/xiaohan2/Software/CREfinder/Demo/02-output/Integration &&

echo "end at `date`" &&
echo "Still water run deep" 1>&2 &&
echo "Still water run deep" > demo.sign
