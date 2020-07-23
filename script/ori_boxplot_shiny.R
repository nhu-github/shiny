#' @import ggpubr
#' 
#' Drawing box plot
#'
#' @param input The input can be three type. One is the mutation record by 0 and 1. The other is origimed var format.The three is clinical data.
#' @param x  x axix type to be show.eg:Mutant gene
#' @param y  y axix type to be show.The y is continuous variable.eg: TMB
#' @param xtype  x from the mutant gene or clinical.eg:"gene" OR "TMB
#' @param cutoff if the x is continuous variable,the cutoff is necessary.
#' @param Output Whether to output pdf file.Default is True.
#'
#' @return A dataframe for ggboxplot
#' @examples
#' ori_boxplot(info1,"Stge分期","TMB",xtype="clinical")
#' ori_boxplot(info1,"年龄","TMB",xtype="clinical")
#' ori_boxplot("TP53","TMB",Mut0,,xtype="gene")
#'
#' @export


## v1 ,input one file 
# ori_boxplot <- function(x, y, xtype, input, ftype="boxplot", cutoff=65, outdir=NULL){
#   require(ggpubr)
#   if(!is.null(outdir)){
#     outdir <- paste0("./",outdir,"/")
#     if(!dir.exists(outdir)){
#       dir.create(outdir)
#     }
#     
#   }
#   if(xtype=="gene"){
#     D  <- dim(table(as.factor(input[1,])))
#     if(D==2){
#       datfra <- data.frame(t(input[c(x,y),]))
#       datfra$type[datfra[[x]]==1] = "MT"
#       datfra$type[datfra[[x]]==0] = "WT"
#     }else{
#       datfra0 <-  unique(input[, c("ORDER_ID", "GENE", y)])
#       datfra <- unique(input[, c("ORDER_ID", y)])
#       MTsample <-unique(datfra0[["ORDER_ID"]][datfra0[["GENE"]]==x])
#       WTsample <- unique(datfra0[["ORDER_ID"]][!datfra0[["ORDER_ID"]] %in% MTsample])
#       datfra$type[datfra$ORDER_ID %in% MTsample ] = "MT"
#       datfra$type[datfra$ORDER_ID %in% WTsample]  = "WT"
#       datfra[[y]] <- as.numeric(datfra[[y]])
#       datfra <-datfra[!is.na(datfra[[y]]),]
#     }
#   }
#   if(xtype=="clinical"){
#     datfra <- input[,c(x,y)]
#     nlevel <- length(levels(as.factor(datfra[[x]])))
#     if(nlevel > 10){
#       datfra <- data.frame(t(apply(datfra, 1, as.numeric)))
#       datfra <- na.omit(datfra)
#       colnames(datfra) <- c(x,y)
#       datfra[[y]] <- as.numeric(datfra[[y]])
#       datfra$type <- ifelse(datfra[[x]] >= cutoff,paste0(x,"_high"),paste0(x,"_low"))
#     }else{
#       if(any(toupper(datfra[[x]]) %in% c("UK","UNKNOW","NULL"))){
#         stop('Please check the record in clinical data')
#       }
#       datfra[[y]] <- as.numeric(datfra[[y]])
#       datfra$type <- as.factor(datfra[[x]])
#     }
#   }
#   
#   if(ftype=="boxplot"){
#     p <- ggboxplot(datfra, x = "type", y = y,
#                    color = "type", palette = "jco",add = "jitter",title = x)+ 
#       stat_compare_means()+
#       #stat_compare_means(label.y =min(datfra[[y]])*1.1) +# Add pairwise comparisons p-value
#       #stat_compare_means(label.y = 0.7) +
#       theme(legend.position='none')+
#       theme(title =  element_text(colour = 'black', angle = 0,size = 18))+
#       theme(plot.title = element_text(hjust = 0.5))
#   }
#   if(ftype=="violin"){
#    p <- ggplot(datfra, aes(x=type, y=datfra[[y]],fill=type))+geom_violin(width=0.8) +
#       geom_boxplot(width=0.05, color="black", fill='white')+
#       ylab(y)+
#      ggtitle(x)+
#       stat_compare_means()+
#       #stat_compare_means(label.y =min(datfra[[y]])*1.1) +# Add pairwise comparisons p-value
#       #stat_compare_means(label.y = 0.7) +
#       theme(legend.position='none')+
#       theme(title =  element_text(colour = 'black', angle = 0,size = 18))+
#       theme(plot.title = element_text(hjust = 0.5))
# 
#     # p <- ggviolin(datfra, x = "type", y = y,
#     #               color = "type", palette = "jco",add = "jitter")+ 
#     #   stat_compare_means()+
#     #   #stat_compare_means(label.y =min(datfra[[y]])*1.1) +# Add pairwise comparisons p-value
#     #   #stat_compare_means(label.y = 0.7) +
#     #   theme(legend.position='none')+
#     #   theme(title =  element_text(colour = 'black', angle = 0,size = 18))+
#     #   theme(plot.title = element_text(hjust = 0.5))
#     # 
#   }
#   
#   #ggsave(paste0(outdir,x,"_vs_",y,".pdf"))
#   return(p)
# }
# 
# 
# 
# # here is change order
# # # example
# # library("openxlsx")
# # library("data.table")
# # 
# # dat <- read.table("RNAseq_Songhui_batch_removal_and_best_response_CR_stat_TMBcount_snvindel_AFF_TMB_batch.txt")
# # Mut0  <- read.xlsx("E:/song/project/20191028ESCC/datbackup.xlsx",sheet = 2)
# # info0 <- read.xlsx("E:/song/project/20191028ESCC/datbackup.xlsx",sheet = 1)
# # 
# # input <- dat[grepl("mut",rownames(dat)),]
# # rownames(input) <- gsub(".*_","",rownames(input))
# # input <- rbind(input,dat["TMB",])
# # input["TMB",]
# # 
# # t <- ori_boxplot("TP53","TMB",xtype="gene",input)
# # t <- ori_boxplot("TP53","TMB",xtype="gene",Mut0)
# # 
# # t <- ori_boxplot("Stge分期","TMB",xtype="clinical",info1)
# # t <- ori_boxplot("age","TMB",xtype="clinical",info1)
# # ori_boxplot("TP53","TMB",xtype="gene",input=datatmp1,ftype="violin", cutoff=65, outdir=NULL)



# v2 ,input 2 files,(mut and clinical)

ori_boxplot <- function(x, y, xtype, input, cli, ftype="boxplot", cutoff=65, outdir=NULL){
  require(ggpubr)
  if(!is.null(outdir)){
    outdir <- paste0("./",outdir,"/")
    if(!dir.exists(outdir)){
      dir.create(outdir)
    }
  }
  if(xtype=="gene"){
    D  <- dim(table(as.factor(input[1,])))
    if(D==2){
      datfra <- data.frame(t(input[c(x,y),]))
      datfra$type[datfra[[x]]==1] = "MT"
      datfra$type[datfra[[x]]==0] = "WT"
      datfra$type <- factor(datfra$type,levels = c("MT","WT"))
    }else{
      datfra0 <-  unique(input[, c("ORDER_ID", "GENE")])
      datfra <- unique(cli[, c("ORDER_ID", y)])
      MTsample <-unique(datfra0[["ORDER_ID"]][datfra0[["GENE"]]==x])
      WTsample <- unique(datfra0[["ORDER_ID"]][!datfra0[["ORDER_ID"]] %in% MTsample])
      datfra$type[datfra$ORDER_ID %in% MTsample ] = "MT"
      datfra$type[datfra$ORDER_ID %in% WTsample]  = "WT"
      datfra$type <- factor(datfra$type,levels = c("MT","WT"))
      datfra[[y]] <- as.numeric(datfra[[y]])
      datfra <-datfra[!is.na(datfra[[y]]),]
    }
  }
  if(xtype=="clinical"){
    datfra <- input[,c(x,y)]
    nlevel <- length(levels(as.factor(datfra[[x]])))
    if(nlevel > 10){
      datfra <- data.frame(t(apply(datfra, 1, as.numeric)))
      datfra <- na.omit(datfra)
      colnames(datfra) <- c(x,y)
      datfra[[y]] <- as.numeric(datfra[[y]])
      datfra$type <- ifelse(datfra[[x]] >= cutoff,paste0(x,"_high"),paste0(x,"_low"))
    }else{
      if(any(toupper(datfra[[x]]) %in% c("UK","UNKNOW","NULL"))){
        stop('Please check the record in clinical data')
      }
      datfra[[y]] <- as.numeric(datfra[[y]])
      datfra$type <- as.factor(datfra[[x]])
    }
  }
  
  if(ftype=="boxplot"){
    p <- ggboxplot(datfra, x = "type", y = y,
                   color = "type", palette = "jco",add = "jitter",title = x)+ 
      stat_compare_means()+
      #stat_compare_means(label.y =min(datfra[[y]])*1.1) +# Add pairwise comparisons p-value
      #stat_compare_means(label.y = 0.7) +
      theme(legend.position='none')+
      theme(title =  element_text(colour = 'black', angle = 0,size = 18))+
      theme(plot.title = element_text(hjust = 0.5))
  }
  if(ftype=="violin"){
    p <- ggplot(datfra, aes(x=type, y=datfra[[y]],fill=type))+geom_violin(width=0.8) +
      geom_boxplot(width=0.05, color="black", fill='white')+
      ylab(y)+
      ggtitle(x)+
      stat_compare_means()+
      #stat_compare_means(label.y =min(datfra[[y]])*1.1) +# Add pairwise comparisons p-value
      #stat_compare_means(label.y = 0.7) +
      theme(legend.position='none')+
      theme(title =  element_text(colour = 'black', angle = 0,size = 18))+
      theme(plot.title = element_text(hjust = 0.5))
    
    # p <- ggviolin(datfra, x = "type", y = y,
    #               color = "type", palette = "jco",add = "jitter")+ 
    #   stat_compare_means()+
    #   #stat_compare_means(label.y =min(datfra[[y]])*1.1) +# Add pairwise comparisons p-value
    #   #stat_compare_means(label.y = 0.7) +
    #   theme(legend.position='none')+
    #   theme(title =  element_text(colour = 'black', angle = 0,size = 18))+
    #   theme(plot.title = element_text(hjust = 0.5))
    # 
  }
  
  #ggsave(paste0(outdir,x,"_vs_",y,".pdf"))
  return(p)
}



