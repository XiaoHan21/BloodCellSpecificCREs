#!/usr/bin/sh
## This is a demo for Software SpecificCREsIdentification
## Written by XiaoHan

echo "start at `date`" &&
<<'COMMENTS'
/jdfsbjcas1/ST_BJ/P21H28400N0232/xiaohan2/Software/miniconda/envs/R4.2/bin/Rscript CREfinder.r \
--Mode RNA \
--ExpMatrix /jdfsbjcas1/ST_BJ/P21H28400N0232/xiaohan2/Software/CREfinder/Demo/01-input/RNA/i1.RNA_TPMat.txt \ 
--SampleMeta /jdfsbjcas1/ST_BJ/P21H28400N0232/xiaohan2/Software/CREfinder/Demo/01-input/RNA/i1.RNA_SampleMeta.txt \
--Features /jdfsbjcas1/ST_BJ/P21H28400N0232/xiaohan2/Software/CREfinder/Demo/01-input/RNA/i1.RNA_DEGs.txt \
--TargetCelltype /jdfsbjcas1/ST_BJ/P21H28400N0232/xiaohan2/Software/CREfinder/Demo/01-input/RNA/i1.RNA_ConcernCelltype.txt \
--Outdir /jdfsbjcas1/ST_BJ/P21H28400N0232/xiaohan2/Software/CREfinder/Demo/02-output/RNA &&

COMMENTS

/jdfsbjcas1/ST_BJ/P21H28400N0232/xiaohan2/Software/miniconda/envs/R4.2/bin/Rscript CREfinder.r \
--Mode ATAC \
--ExpMatrix /jdfsbjcas1/ST_BJ/P21H28400N0232/xiaohan2/Software/CREfinder/Demo/01-input/ATAC/i1.ATAC_CPMat.txt \
--SampleMeta /jdfsbjcas1/ST_BJ/P21H28400N0232/xiaohan2/Software/CREfinder/Demo/01-input/ATAC/i1.ATAC_SampleMeta.txt \
--Features /jdfsbjcas1/ST_BJ/P21H28400N0232/xiaohan2/Software/CREfinder/Demo/01-input/ATAC/i1.ATAC_DORs.txt \
--TargetCelltype /jdfsbjcas1/ST_BJ/P21H28400N0232/xiaohan2/Software/CREfinder/Demo/01-input/ATAC/i1.ATAC_ConcernCelltype.txt \
--Outdir /jdfsbjcas1/ST_BJ/P21H28400N0232/xiaohan2/Software/CREfinder/Demo/02-output/ATAC &&

echo "end at `date`" &&
echo "Still water run deep" 1>&2 &&
echo "Still water run deep" > demo.sign
