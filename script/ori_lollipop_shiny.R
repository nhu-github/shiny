#' @include lollipop.r
#' @import maftools
NULL
#' Drawing lollipop plot
#'
#' @param gene  The gene to draw lollipop plot
#' @param mut Mutation file
#' @param refSeqID RefSeq transcript identifier for \code{gene} if known.Notice the string split by ";"
#' @param outdir The output directory name
#'
#' @return 
#' @examples
#' ori_lollipop("TP53",Mut1,"refSeqID","lollipop")
#'
#'
#' @export

ori_lollipop <- function(gene,mut,refSeqID,outdir=NULL){
  if(!is.null(outdir)){
    outdir <- paste0("./",outdir,"/")
    if(!dir.exists(outdir)){
      dir.create(outdir)
    }
  }
  gff <- readRDS("./data/gff.rds")
  #maf <- maf[which(maf$REPORT_OR_VUS=="report"),]
  maf <-  mut[which(mut$GENE == gene), ]
  refID <- sapply(strsplit(maf[[refSeqID]],";"),"[",1) # notice the ";"
  refID <-  names(sort(table(refID),decreasing = T))[1]
  if(refID %in% gff$refseq.ID){
    #pdf(paste0(outdir,"lollipop_",gene,".pdf"))
    lollipop(maf = maf, gene = gene, gff=gff, AACol = 'AA_CHANGE', refSeqID = refID)
    #dev.off()
  }
}

# rm(list = ls())
# 
# library("data.table")
# library("maftools")
# source("E:/song/project/20191106biomarker/data/wyldata/script/lollipop.r")
# dir <- "E:/song/project/20191106biomarker/data/wyldata/script/gff.RData"
# load(dir)
# Mut1 <- read.xlsx("E:/song/project/20191106biomarker/data/wyldata/longhuamut.xlsx")
# # setnames(Mut1,"参考序???","refSeqID")
# # Mut1$refSeqID <- gsub("??",";",Mut1$refSeqID )




