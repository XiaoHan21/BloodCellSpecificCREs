#!/usr/bin/sh
## This is a demo for Software MockReadsGenerator
## Written by XiaoHan
## 2023-09-20

echo "start at `date`" &&
/jdfsbjcas1/ST_BJ/P21H28400N0232/xiaohan2/Software/miniconda/envs/R4.2/bin/Rscript step1.MockReadsGenerator.R \
--Mode BaseRandom \
--Templates ./Demo/01-input/BaseRandom/i1.refseq.txt \
--SiteNums 300 \
--ReadLength 100 \
--NegativeRate 0.1 \
--BarcodeLens 8 \
--Outdir Demo/02-output/BaseRandom &&

/jdfsbjcas1/ST_BJ/P21H28400N0232/xiaohan2/Software/miniconda/envs/R4.2/bin/Rscript step1.MockReadsGenerator.R \
--Mode TemplateRandom \
--Templates ./Demo/01-input/TemplateRandom/i1.RefSeqs.xlsx \
--SiteNums 500 \
--ReadLength 150 \
--Outdir Demo/02-output/TemplateRandom &&

bash step2.MakeFqFormat.sh Demo/02-output/BaseRandom/o1.BaseRandom_demo.fq Demo/02-output/BaseRandom/o1.BaseRandom.fq &&
bash step2.MakeFqFormat.sh Demo/02-output/TemplateRandom/o1.TemplateRandom_demo.fq Demo/02-output/TemplateRandom/o1.TemplateRandom.fq &&
echo "end at `date`" &&
echo "Still water run deep" 1>&2 &&
echo "Still water run deep" > demo.sign
