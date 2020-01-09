#' @import ggpubr
NULL
#' Drawing cor plot
#'
#' @param data The rownames as sample name, the column as feature name,record with continuous values
#' @param x The colnames of data as x axix,can be one or more
#' @param y The colnames of data as y axix,only one 
#' @param outdir the outdir name
#' @param method a character string indicating which correlation coefficient (or covariance) is to be computed. One of "pearson" (default), "kendall", or "spearman".
#'
#' @return A ggplot2 object or pdf
#' @examples
#' ori_scattercor_plot(datacor,x= c("Monocytes","Neutrophils"),y="TMB")
#' ori_scattercor_plot(datacor,x= "Monocytes",y="TMB")
#'
#' @export

ori_scattercor_plot <- function(x,y,data,method="pearson",outdir=NULL){
  require(ggpubr)
  if(!is.null(outdir)){
    outdir <- paste0("./",outdir,"/")
    if(!dir.exists(outdir)){
      dir.create(outdir)
    }
    
  }
  L <- length(c(x,y))
  if(L==2){
    p <-ggscatter(data, x = x, y = y, size = 0.3,palette = "jco",
                  add = "reg.line", conf.int = TRUE)+
      stat_cor(method = method)
  }
  if(L>2){
    Co <- c(y,x)
    data <- data[,Co]
    data <- tidyr::gather(data, key = 'Variants', value = 'value',-y)
    data[[1]] <- as.numeric(data[[1]])
    data$value <- as.numeric(data$value)
    p <- ggscatter(data, x = y, y = "value", size = 0.3,
                   color = 'Variants', palette = "jco",
                   facet.by = 'Variants', #scales = "free_x",
                   add = "reg.line", conf.int = TRUE) +
      stat_cor(aes(color = Variants), method = method)
  }
  #ggsave(plot = p, filename = paste0(outdir,x,"_vs_",y,".pdf"), height = 9, width = 6)
  return(p)
}


# dat <- read.table("RNAseq_Songhui_batch_removal_and_best_response_CR_stat_TMBcount_snvindel_AFF_TMB_batch.txt")
# dat1 <- dat[1:27864,]
# write.table(dat1,"expression.txt",sep="\t")
# source('CIBERSORT.R')
# exprdir <- "./expression.txt"
# LM22dir <- "./LM22.txt"
# ciber_result <- CIBERSORT(exprdir,LM22dir, perm = 200, absolute = F)
# ciber_result1 <- as.data.frame(t(ciber_result))
# ciber_result1 <- ciber_result1[1:(nrow(ciber_result1)-3),]
# TMB <- as.data.frame(t(dat["TMB",]))
# TMB <- TMB[match(rownames(ciber_result1),rownames(TMB))]
# datacor <- cbind(ciber_result1,TMB)
# ori_scattercor_plot(datacor,x= c("Monocytes","Neutrophils"),y="TMB")
# ori_scattercor_plot(datacor,x= "Monocytes",y="TMB")