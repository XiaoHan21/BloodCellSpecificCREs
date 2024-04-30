########################################################################################
# Author: xiaohan2@genomics.cn
# Date: 2023-09-20 12:00:00
# LastEditors: xiaohan2@genomics.cn
# LastEditTime: 2023-09-20 12:00:00
# Description: 生成抓取N8 barcode的参数测试集
# Mock dataset type I: 
#                       目的：滑动生成覆盖某个序列区间的mock reads（序列区间碱基随机生成）
#                       输入：模板参考序列、序列区间长度、阴性对照reads率 
#                       输出：一组覆盖指定区间序列的mock reads，指定区间序列随机
# Mock dataset type II: 
#                       目的：滑动生成覆盖某个序列区间的mock reads（序列区间碱基通过指定）
#                       输入：一组模板参考序列、序列区间长度
#                       输出：一组覆盖指定区间序列的mock reads，指定区间序列从参考模板选取
# Test dataset1: 从barcode侧翼位点滑动的91个位点的100bp reads,
#                其中每个相同位点300条reads，270条含有随机barcode,30条不含随机barcode
# Test dataset2: 从39条参考序列区间中随机挑选一条，滑动149个位点的150bp reads，
#                其中每个相同位点500条reads
# Copyright 2023 by xiaohan2 All Rights Reserved.
########################################################################################

### >>> 0. 参数传递
#-------------------------------------------------------------------------------
library(optparse)
library(openxlsx)
library(stringr)
library(dplyr)

# 描述参数
option_list <- list(
  make_option(c("--Mode"), type = "character", default = FALSE,
              action = "store", help = "Mode: BaseRandom or TemplateRandom"
  ),
  make_option(c("--Templates"), type = "character", default = FALSE,
              action = "store", help = "MotherSeq templates file path"
  ),
  make_option(c("--SiteNums"), type = "integer", default = FALSE,
              action = "store", help = "The number of single site mock reads"
  ),
  make_option(c("--ReadLength"), type = "integer", default = FALSE,
              action = "store", help = "The length of mock reads"
  ),
  make_option(c("--NegativeRate"), type = "double", default = FALSE,
              action = "store", help = "ONLY for BaseRandom mode: the number of single site negative mock reads"
  ),
  make_option(c("--BarcodeLens"), type = "integer", default = FALSE,
              action = "store", help = "ONLY for BaseRandom mode: the length of barcodes"
  ),
  make_option(c("--Outdir"), type = "character", default = FALSE,
              action = "store", help = "Result directory"
  ),
  make_option(c("--Try"), type = "logical", default = FALSE,
              action = "store_TRUE", help = "This is demo and have a try!"
  )
)

# 解析参数
args = parse_args(OptionParser(option_list = option_list, usage = "This script is a test for arguments!"))

### >>> Main <<<
#-------------------------------------------------------------------------------
### Part 1: BaseRandom
if(args$Mode=="BaseRandom"){
  if (args$Try) { # 程序测试使用
    outdir <- "Demo/02-output/BaseRandom"
    if(!dir.exists(outdir)){dir.create(outdir,recursive = T)}
    # 模版序列需替换位置用N表示
    mother_seq <- "cctggagacctccgcgccccgcaacctccccttctacgagcggctcggcttcaccgtcaccgccgacgtcgaggtgcccgaaggaccgcgcacctggtgcatgacccgCAAGCCCGGTGCCTGATGCAGGCATATCAATAAGCGGAGGANNNNNNNNCGATATCTCGAGGGTACCTTTAAGACCAATGACTTACAAGGCAGCTGTAGATCTTAGCCACTTTTTAAAAGAAAAGGGGGGACTGGAAGGGCTAATTCACTCCCAACGAAGATAAGATCTGCTTTTTGCTTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGG"
    # mother_seq <- read.table("Demo/01-input/BaseRandom/i1.refseq.txt",header=F)$V1
    n_barcode_len <- 8 # Barcode长度
    n_site_reads <- 300 # 单个位点总reads数
    neg_rate <- 0.1 # 单个位点阴性对照数
    read_lens <- 100 # 单个reads长度
    barcode.tag <- paste(rep("N",n_barcode_len),collapse="") # 模版序列中的N
  } else {
    outdir <- args$Outdir
    if(!dir.exists(outdir)){dir.create(outdir,recursive = T)}
    # 模版序列需替换位置用N表示
    mother_seq <- read.table(args$Templates,header=F)$V1 # 模板序列
    n_barcode_len <- args$BarcodeLens # Barcode长度
    n_site_reads <- args$SiteNums # 单个位点总reads数
    neg_rate <- args$NegativeRate # 单个位点阴性对照数
    read_lens <- args$ReadLength # 单个reads长度
    barcode.tag <- paste(rep("N",n_barcode_len),collapse="") # 模版序列中的N
  }
  ### >>> 1. 确定100bp目标序列(分带/不带barcode)
  #-------------------------------------------------------------------------------
  # 确定Barcode起始位点与终末位点
  N_index_begin <- gregexpr("N",mother_seq)[[1]][1]
  N_index_end <- gregexpr("N",mother_seq)[[1]][n_barcode_len]
  
  # 滑动起始终止位点
  start_site_index <- N_index_end + 1 - (read_lens - 1) 
  end_site_index <- N_index_begin - 1
  
  ### >>> 2. 生成随机字串
  #-------------------------------------------------------------------------------
  N_random_barcode <- function(N_len,N_reads,sed=666){
    set.seed(sed)
    for (i in 1:N_len) {
      if (i == 1) {
        out <- sample(c("A","G","C","T"),N_reads,replace = T)
      } else{
        out <- paste0(out, sample(c("A","G","C","T"),N_reads,replace = T))
      }
    }
    return(out)
  }
  
  ### >>> 3. 位点字符串提取
  #-------------------------------------------------------------------------------
  site_seq <- function(mother_str,start_site,n_len_reads,n_len_barcode,n_site_reads,neg_rate){
    ### 提取目标区段序列
    pos_seq <- substr(mother_str,start_site,start_site + n_len_reads - 1)
    neg_seq <- gsub("N","",substr(mother_str,start_site - n_len_barcode,start_site + n_len_reads - 1)) # 阴性对照reads与阳性reads终止位点相同，起始位点平移一个barcode长度
    ### Barcode在序列中的位置
    barcode_start_site <- unlist(gregexpr(barcode.tag,pos_seq,fixed = T))
    barcode_end_site <- barcode_start_site + n_len_barcode - 1
    ### 计算两组reads数目
    pos_num <- n_site_reads * (1 - neg_rate)
    neg_num <- n_site_reads * neg_rate
    ### 根据pos reads数产生相应数目的随机barcode
    sub_pos_randomBar <- N_random_barcode(N_len = n_len_barcode,N_reads = pos_num,sed = start_site)
    ### 重复目标区段
    pos_seq_n <- rep(pos_seq,pos_num)
    neg_seq_n <- rep(neg_seq,neg_num)
    bar_start_site_n <- rep(barcode_start_site,pos_num)
    bar_end_site_n <- rep(barcode_end_site,pos_num)
    ### 替换上随机barcode
    for (i in 1:pos_num) {
      pos_seq_n[i] <- gsub(barcode.tag,sub_pos_randomBar[i],pos_seq_n[i])
    }
    
    ### 结果整合
    out.list <- list()
    barcode.seq <- c(sub_pos_randomBar,rep(NA,neg_num))
    barcode.start.index <- c(bar_start_site_n,rep(NA,neg_num))
    barcode.end.index <- c(bar_end_site_n,rep(NA,neg_num))
    total_seq <- c(pos_seq_n,neg_seq_n)
    out.list <- list(total_seq,barcode.seq,barcode.start.index,barcode.end.index)
    return(out.list)
  }
  
  for (i in start_site_index:end_site_index) {
    tmp <- site_seq(mother_str = mother_seq,
                    n_len_reads = read_lens,
                    n_site_reads = n_site_reads,
                    start_site = i,
                    n_len_barcode = n_barcode_len,
                    neg_rate = neg_rate)
    
    seq <- tmp[[1]]
    barcode.seq <- tmp[[2]]
    barcode.start.index <- tmp[[3]]
    barcode.end.index <- tmp[[4]]
    
    if(i == start_site_index){
      out.seq <- seq
      out.barcode.seq <- barcode.seq
      out.barcode.start.site <- barcode.start.index
      out.barcode.end.site <- barcode.end.index
    } else{
      out.seq <- c(out.seq,seq)
      out.barcode.seq <- c(out.barcode.seq,barcode.seq)
      out.barcode.start.site <- c(out.barcode.start.site,barcode.start.index)
      out.barcode.end.site <- c(out.barcode.end.site,barcode.end.index)
    }
  }
  
  ### >>> 4. 输出为fq文件
  #-------------------------------------------------------------------------------
  total_reads_num <- (end_site_index - start_site_index + 1) * n_site_reads
  reads_id <- str_pad(1:total_reads_num,width = 7,side = "left",pad = "0")
  reads_id <- paste0("@RI:",reads_id)
  barcode_id <- paste0("BC:",gsub(">","",out.barcode.seq))
  start_id <- paste0("SS:",out.barcode.start.site)
  end_id <- paste0("ES:",out.barcode.end.site)
  
  ## 整理fastq文件
  info.merge <- paste(reads_id,barcode_id,start_id,end_id,sep = "|||")
  
  skeleton.fq <- data.frame(info=info.merge,seq=out.seq)
  skeleton.fq <- skeleton.fq %>% mutate(strand="+",score=paste0(rep("F",read_lens),collapse = ""))
  
  ## 整理metadata文件
  metadata <- data.frame(seqname=info.merge, seq=out.seq, barcode=out.barcode.seq,
                         barcode_start_site=out.barcode.start.site, barcode_end_site=out.barcode.end.site)
  metadata <- metadata %>% mutate(is_empty=case_when(grepl("NA",metadata$seqname)~"WithoutBarcode",
                                                     TRUE~"WithBarcode"))
  
  ## 输出序列
  write.table(skeleton.fq,paste0(outdir,"/o1.BaseRandom_demo.fq"),sep = "\t",quote = F,col.names = F,row.names = F)
  write.xlsx(metadata,paste0(outdir,"/o2.BaseRandom_metadata.xlsx"),colNames=T,rowNames=F,keepNA=T)

### Part 2: TemplateRandom
}else if (args$Mode=="TemplateRandom"){
  if (args$Try) { # 程序测试使用
    ### >>> 1. 配置文件
    ###-----------------------------------------------------------------------------
    outdir <- "Demo/02-output/TemplateRandom"
    if(!dir.exists(outdir)){dir.create(outdir,recursive = T)}
    source.dt <- read.xlsx("Demo/01-input/TemplateRandom/i1.RefSeqs.xlsx",sheet = 1,startRow = 1)
    n_site_reads <- 500 # 单个位点总reads数
    read_lens <- 150 # 单个reads长度
    # 滑动起始终止位点
    start_site_index <- 1
    end_site_index <- 149
  } else{
    ### >>> 1. 配置文件
    ###-----------------------------------------------------------------------------
    outdir <- args$Outdir
    if(!dir.exists(outdir)){dir.create(outdir,recursive = T)}
    source.dt <- read.xlsx(args$Templates,sheet = 1,startRow = 1)
    n_site_reads <- args$SiteNums
    read_lens <- args$ReadLength
    start_site_index <- 1
    end_site_index <- read_lens - 1
  }
  
  ### >>> 2. 位点字符串提取函数
  #-------------------------------------------------------------------------------
  ### 对于每个位点，500条reads，每条reads长度150bp，随机从39种参考序列中抽取
  site_seq_random <- function(meta,start_site,n_len_reads,n_site_reads,seeds=666){
    ### 原始信息提取
    pos_num <- nrow(meta)
    pos_str <- meta[,2]
    pos_seqname <- meta[,1]
    
    set.seed(seed = seeds)
    seed.num <- sample(1:pos_num,n_site_reads,replace = T)
    
    ### 提取目标区段序列
    pos_seq <- substr(pos_str[seed.num],start_site,start_site + n_len_reads - 1)
    pos_seqname_record <- pos_seqname[seed.num]
    
    ### 结果整合
    total_seq <- pos_seq
    total_seqname <- pos_seqname_record
    site_start <- rep(start_site,n_site_reads)
    site_end <- rep(start_site + n_len_reads - 1,n_site_reads)
    out.list <- list(total_seqname,total_seq,site_start,site_end)
    return(out.list)
  }
  
  ### >>> 3. 滑动生成数据
  #-------------------------------------------------------------------------------
  for (i in start_site_index:end_site_index) {
    tmp <- site_seq_random(meta = source.dt[,1:2],
                           start_site = i,
                           n_len_reads = read_lens,
                           n_site_reads = n_site_reads,
                           seeds = i)
    
    seqname <- tmp[[1]]
    seq <- tmp[[2]]
    start.index <- tmp[[3]]
    end.index <- tmp[[4]]
    
    if(i == start_site_index){
      out.seqname <- seqname
      out.seq <- seq
      out.start.index <- start.index
      out.end.index <- end.index
    } else{
      out.seqname <- c(out.seqname,seqname)
      out.seq <- c(out.seq,seq)
      out.start.index <- c(out.start.index,start.index)
      out.end.index <- c(out.end.index,end.index)
    }
  }
  
  ### >>> 4. 输出为fq文件
  #-------------------------------------------------------------------------------
  total_reads_num <- (end_site_index - start_site_index + 1) * n_site_reads
  reads_id <- str_pad(1:total_reads_num,width = 7,side = "left",pad = "0")
  reads_id <- paste0("@RI:",reads_id)
  barcode_id <- paste0("BC:",gsub(">","",out.seqname))
  start_id <- paste0("SS:",out.start.index)
  end_id <- paste0("ES:",out.end.index)
  
  info.merge <- paste(reads_id,barcode_id,start_id,end_id,sep = "|||")
  
  skeleton.fq <- data.frame(info=info.merge,seq=out.seq)
  skeleton.fq <- skeleton.fq %>% mutate(strand="+",score=paste0(rep("F",read_lens),collapse = ""))
  
  metadata <- data.frame(seqname=info.merge,seq=out.seq,start_site=out.start.index,end_site=out.end.index)
  metadata <- metadata %>% mutate(is_empty=case_when(grepl("EmptyVector",metadata$seqname)~"empty",
                                                     TRUE~"Non_empty"))
  
  write.table(skeleton.fq,paste0(outdir,"/o1.TemplateRandom_demo.fq"),
              sep = "\t",quote = F,col.names = F,row.names = F)
  write.xlsx(metadata,paste0(outdir,"/o2.TemplateRandom_metadata.xlsx"),colNames=T,rowNames=F,keepNA=T)
}
