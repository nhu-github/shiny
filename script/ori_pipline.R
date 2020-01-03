#' @include CountVariants.R
#' @include CIBERSORT.R
#' @include lollipop.r
#' @include ori_profiling.R 
#' @include ori_boxplot.R   
#' @include ori_cibersort.R 
#' @include ori_mut_signature.R
#' @include ori_scattercor_plot.R
#' @include ori_pheatmap.R
#' @include ori_biomarker_gene.R
#' @include ori_CountVariant.R
#' @include ori_barplot_draw.R 
#' @include ori_barplot.R 
#' @include ori_tools.R 
#' @import openxlsx
#' @import data.table
#' @import ComplexHeatmap
#' @import ggplot2
#' @import deconstructSigs
#' @import BSgenome.Hsapiens.1000genomes.hs37d5
#' @import ggsci
#' @import pheatmap
#' @import RColorBrewer
#' @import maftools
#' 
#' biomarker pipline
#'
#' @param mut Mutation file, there must be one colnames named "ORDER_ID"
#' @param cli Clinical file
#' @param exprdir The directory where stored the expression file. The field separator character is "\t"
#' @param bar Column names of clinical file used for draw barplot on the top of profiling.Continous variable.
#' @param feature Column names of clinical file used for pheatmap on the top of profiling.Discrete variable.
#' @param feature1 Column names of clinical file used for statistics
#' @param prefix output pdf file's prefix
#' @param pathway_genes Pathway dataframe, column as pathway name,rownames as the gene names in the pathway
#' @param n The number of first top mutantdgenes to show.Default is 30.
#' @param nfreq The cutoff of mutation frequency to show.
#' @param cutoff The cutoff of groupby.Default is the median of groupby value.
#'
#' @return output the pdf files
#' @examples
#'  load("Demo.Rdata")
#'  ori_pipline(mut=Mut1,cli=info1,expdir="./expression1.txt",bar="TMB",
#'  feature=c("GENDER","Smoker"),prefix="longhua",pathway_genes=pathway)
#'
#' @export


ori_pipline <- function(mut, 
                        cli,
                        expdir=NULL,
                        bar,
                        feature,
                        feature1,
                        prefix, 
                        pathway_genes=NULL, 
                        n = 30,
                        nfreq=0.03, 
                        cutoff=NULL){
  options(warn=-1)
  
  mut <- openxlsx::read.xlsx(mut, sheet = 1, skipEmptyRows = TRUE,
                                               na.strings = c("NA", "NaN"), detectDates = TRUE)
  cli <- openxlsx::read.xlsx(cli, sheet = 1, skipEmptyRows = TRUE,
                             na.strings = c("NA", "NaN"), detectDates = TRUE)
  if(!is.null(pathway_genes)){pathway_genes<- openxlsx::read.xlsx(pathway_genes)}

  stat <- ori_biomarker_gene(bar,mut,n=5,cutoff=cutoff,outdir="stat")
  statcli <- ori_statcli(feature1,bar,cli,cutoff=cutoff,outdir="stat")
  CountVariant <-ori_CountVariant(mut,bar,cutoff = cutoff,outdir="stat")
  barp <- ori_barplot_draw(mut,bar,cutoff = cutoff,outdir="barplot")
  landscap <- plot_landscape(mut,cli,bar,feature,prefix,cutoff=cutoff,n = n)
  if(!is.null(pathway_genes)){
    landscap_path <- plot_pathway_cli(mut, cli, bar, feature, prefix, pathway_genes, nfreq=nfreq, cutoff=cutoff)
  }
  
  freq <- cal_freq(mut, length(unique(mut[["ORDER_ID"]])))
  gene <- as.character(freq$GENE[freq$MUTANT>5])
  if(length(gene)> 100){
    gene <- as.character(freq$GENE)[1:100]
  }else{
    gene <- as.character(freq$GENE[freq$MUTANT>5])
  }
  boxp <- lapply(gene,ori_boxplot,bar,"gene",mut,outdir="gene_vs_biomarker")
  
  # ori_lollipop
  if("GENOMIC" %in% rownames(mut)){
    mutsig <- ori_mut_signature(mut, groupby= bar, outdir="siganature")
  }
  
  mut1 <- mut[!mut$AA_CHANGE=="-",]
  freq1 <- cal_freq(mut1, length(unique(mut1[["ORDER_ID"]])))
  gene1 <- as.character(freq1$GENE[freq1$MUTANT>3])
  if(length(gene1)> 100){
    gene1 <- as.character(freq1$GENE)[1:100]
  }else{
    gene1 <- as.character(freq1$GENE[freq1$MUTANT>3])
  }
  lolli <- lapply(gene1,ori_lollipop,mut,"refSeqID","lollipop")
  
  # expdir exist 
  if(!is.null(expdir)){
    ciber<- ori_cibersort(expdir)
    bardat <- cli[match(colnames(ciber),cli$ORDER_ID),bar]
    ciber1 <- ciber[,1:(ncol(ciber)-3)]
    datsca <- data.frame(t(rbind(ciber1,bardat)))
    colnames(datsca)[ncol(datsca)] <- bar
    scattercor <- lapply(colnames(datsca)[1:(ncol(datsca)-1)],ori_scattercor_plot,bar,datsca,outdir="immuno")
    # pheatmap
    data <- read.table(expdir,header=T,sep="\t",row.names=1,check.names=F)
    bardat <- cli[match(colnames(data),cli$ORDER_ID),bar]
    phedat <- rbind(bardat,data)
    rownames(phedat)[1] <- bar
    pheatmap <- ori_pheatmap(phedat[1:(n+1),],bar,cutoff=cutoff)
    
  }
  }


# example
# rm(list = ls())
# load=T
# if(load==T){
#   source("E:/song/project/20191106biomarker/data/wyldata/script/ori_profiling.R")
#   source("E:/song/project/20191106biomarker/data/wyldata/script/ori_boxplot.R")
#   source("E:/song/project/20191106biomarker/data/wyldata/script/ori_mut_signature.R")
#   source("E:/song/project/20191106biomarker/data/wyldata/script/ori_cibersort.R")
#   source("E:/song/project/20191106biomarker/data/wyldata/script/ori_scattercor_plot.R")
#   source("E:/song/project/20191106biomarker/data/wyldata/script/ori_pheatmap.R")
#   source("E:/song/project/20191106biomarker/data/wyldata/script/lollipop.r")
#   source("E:/song/project/20191106biomarker/data/wyldata/script/ori_lollipop.r")
#   source("E:/song/project/20191106biomarker/data/wyldata/script/ori_biomarker_gene.R")
#   source("E:/song/project/20191106biomarker/data/wyldata/script/CountVariants.R")
#   source("E:/song/project/20191106biomarker/data/wyldata/script/ori_barplot.R")
#   source("E:/song/project/20191106biomarker/data/wyldata/script/ori_tools.R")
#   source("E:/song/project/20191106biomarker/data/wyldata/script/ori_CountVariant.R")
#   source("E:/song/project/20191106biomarker/data/wyldata/script/ori_barplot_draw.R")
#   source('E:/song/project/20191106biomarker/data/wyldata/script/CIBERSORT.R')
#   source('E:/song/project/20191106biomarker/data/wyldata/script/ori_fisherplot.R')
#   source('E:/song/project/20191106biomarker/data/wyldata/script/ori_statcli.R')
#   source('E:/song/project/20191106biomarker/data/wyldata/script/survival_zsh.R')
# 
# }
# lib <- T
# if(lib==T){
#   library("openxlsx")
#   library("data.table")
#   library("ComplexHeatmap")
#   library("ggplot2")
#   library("deconstructSigs")
#   library("BSgenome.Hsapiens.1000genomes.hs37d5")
#   library("ggsci")
#   library("pheatmap")
#   library("RColorBrewer")
#   library("maftools")
# }
# 
# load("Demo.Rdata")
# save(list = ls(all=T),file= "Demo.Rdata")
# 
# rm(list = ls())
# load("Demo.Rdata")
# setwd("E:/song/project/20191106biomarker/data/wyldata/test2")
# t <-ori_pipline(mut=Mut1,cli=info1,expdir="./expression1.txt",bar="TMB",
#                 feature=c("GENDER","Smoker"),
#                 feature1=c("GENDER","Smoker","Lauren","Hypertension","Surgery","MSS"),
#                 prefix="test",pathway_genes=pathway)







  

