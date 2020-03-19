#' @param x The column of the clinical feature.Discrete variable
#' @param y  The column that  dividing the data into two groups.eg: TMB
#' @param info   clinical data
#' @param cutoff   The cutoff of y
#' @param outdir output directory
#' 
#' @return output a xlsx file
#' @examples
#'x <- c("GENDER","Smoker","Lauren","Hypertension","Surgery","MSS")
#'y <- "TMB"
#'ori_statcli(x,y,info)
#'
#' @export

#====v1======
# ori_statcli <- function(x,y,info,cutoff=NULL,outdir=NULL){
#   # age column is AGE
#   options(digits=3)
#   summaryCli <- function(x,y,dat){
#     xtable <- table(dat[,c(x, y)])
#     xrp <- sweep(xtable, 2, apply(xtable, 2, sum), '/')
#     xrp <- paste0(' (', round(xrp, digits = 3)*100, '%', ')')
#     xr <- xtable
#     xr[,1] <- paste0(xr[,1], xrp[1:nrow(xr)])
#     xr[,2] <- paste0(xr[,2], xrp[1:nrow(xr)+nrow(xr)])
#     if(nrow(xr) > 1){
#       xfisher <- fisher.test(xtable)
#       xr <- cbind(xr, c(round(xfisher$p.value, digits = 3), rep(NA, nrow(xtable)-1)))
#     }else{
#       xr <- cbind(xr, '-')
#     }
#     xr <- as.matrix(xr)
#     xr <- as.data.frame(xr, stringsAsFactors = F)
#     xr[,4] <- paste0('    ', rownames(xr))
#     xr <- xr[,c(4,1:3)]
#     xr <- rbind(c(x, rep(NA, 3)), xr)
#     return(xr)
#   }
#   
#   if(!is.null(outdir)){
#     outdir <- paste0("./",outdir,"/")
#     if(!dir.exists(outdir)){
#       dir.create(outdir)
#     }
#   }
#   
#   if("AGE" %in% x){
#     x <- x[!x %in% "AGE"]
#     mx <- info[,c("AGE",y),drop=F]
#     if(is.null(cutoff)){
#       mx$type <- ifelse(mx[[y]] > median(mx[[y]],na.rm = T),"High","Low") 
#     }else{
#       mx$type <- ifelse(mx[[y]] > cutoff,"High","Low") 
#     }
#     cgc.age <- aggregate(AGE~type, mx, median)
#     cgc.age <- paste(round(as.numeric(cgc.age[[2]]),digits = 3), apply(aggregate(AGE~type, mx, range), 1, function(x){paste0('(', x[2], '-', x[3], ')')}), sep = ' ')
#     cgc.age <- c('Age', cgc.age, wilcox.test(AGE~type, mx)$p.value)
#     
#     info <- cbind(info,"type"=mx$type)
#     cgc.table <- lapply(x, summaryCli,"type",info)
#     cgc.table <- do.call(rbind, cgc.table)
#     
#     cgc.table <- rbind(cgc.age, cgc.table)
#     cgc.table$V3 <- round(as.numeric(cgc.table$V3),digits = 3)
#     colnames(cgc.table)[1] <- "characteristics"
#     colnames(cgc.table)[4] <- "pvalue"
#     #openxlsx::write.xlsx(cgc.table, file = paste0(outdir, 'characteristics.xlsx'), row.names = F)
#     return(cgc.table)
#   }else{
#     mx <- info[,y,drop=F]
#     if(is.null(cutoff)){
#       mx$type <- ifelse(mx[[y]] > median(mx[[y]],na.rm = T),"High","Low") 
#     }else{
#       mx$type <- ifelse(mx[[y]] > cutoff,"High","Low") 
#     }
#     info <- cbind(info,"type"=mx$type)
#     cgc.table <- lapply(x, summaryCli,"type",info)
#     cgc.table <- do.call(rbind, cgc.table)
#     
#     #cgc.table <- rbind(cgc.age, cgc.table)
#     cgc.table$V3 <- round(as.numeric(cgc.table$V3),digits = 3)
#     colnames(cgc.table)[1] <- "characteristics"
#     colnames(cgc.table)[4] <- "pvalue"
#     #openxlsx::write.xlsx(cgc.table, file = paste0(outdir, 'characteristics.xlsx'), row.names = F)
#     return(cgc.table)
# 
#   }
#  
# }

#====v2====
ori_statcli <- function(x,y,info,ytype ="continuous",cutoff=NULL,outdir=NULL){
  # age column is AGE
  options(digits=3)
  summaryCli <- function(x,y,dat){
    xtable <- table(dat[,c(x, y)])
    xrp <- sweep(xtable, 2, apply(xtable, 2, sum), '/')
    xrp <- paste0(' (', round(xrp, digits = 3)*100, '%', ')')
    xr <- xtable
    xr[,1] <- paste0(xr[,1], xrp[1:nrow(xr)])
    xr[,2] <- paste0(xr[,2], xrp[1:nrow(xr)+nrow(xr)])
    if(nrow(xr) > 1){
      xfisher <- fisher.test(xtable)
      xr <- cbind(xr, c(round(xfisher$p.value, digits = 3), rep(NA, nrow(xtable)-1)))
    }else{
      xr <- cbind(xr, '-')
    }
    xr <- as.matrix(xr)
    xr <- as.data.frame(xr, stringsAsFactors = F)
    xr[,4] <- paste0('    ', rownames(xr))
    xr <- xr[,c(4,1:3)]
    xr <- rbind(c(x, rep(NA, 3)), xr)
    return(xr)
  }
  
  if(!is.null(outdir)){
    outdir <- paste0("./",outdir,"/")
    if(!dir.exists(outdir)){
      dir.create(outdir)
    }
  }
  
  if("AGE" %in% x){
    x <- x[!x %in% "AGE"]
    mx <- info[,c("AGE",y),drop=F]
    if(ytype=="continuous"){
      if(is.null(cutoff)){
        mx$type <- ifelse(mx[[y]] > median(mx[[y]],na.rm = T),"High","Low") 
      }else{
        mx$type <- ifelse(mx[[y]] > cutoff,"High","Low") 
      }
    }else{
      mx$type <- as.character(mx[[y]])
    }
    #mx[["AGE"]] <- gsub("岁","",mx[["AGE"]])
    mx[["AGE"]] <- as.numeric(mx[["AGE"]])
    mx <- data.frame(mx)
    cgc.age <- aggregate(AGE~type, mx, median)
    cgc.age <- paste(round(as.numeric(cgc.age[[2]]),digits = 3), apply(aggregate(AGE~type, mx, range), 1, function(x){paste0('(', x[2], '-', x[3], ')')}), sep = ' ')
    cgc.age <- c('Age', cgc.age, wilcox.test(AGE~type, mx)$p.value)
    info <- cbind(info,"type"=mx$type)
    cgc.table <- lapply(x, summaryCli,"type",info)
    cgc.table <- do.call(rbind, cgc.table)
    cgc.table <- rbind(cgc.age, cgc.table)
    cgc.table$V3 <- round(as.numeric(cgc.table$V3),digits = 3)
    colnames(cgc.table)[1] <- "characteristics"
    colnames(cgc.table)[4] <- "pvalue"
    #openxlsx::write.xlsx(cgc.table, file = paste0(outdir, 'characteristics.xlsx'), row.names = F)
    return(cgc.table)
  }else{
    mx <- info[,y,drop=F]
    if(ytype=="continuous"){
      if(is.null(cutoff)){
        mx$type <- ifelse(mx[[y]] > median(mx[[y]],na.rm = T),"High","Low") 
      }else{
        mx$type <- ifelse(mx[[y]] > cutoff,"High","Low") 
      }
    }else{
      mx$type <- as.character(mx[[y]])
    }
    info <- cbind(info,"type"=mx$type)
    cgc.table <- lapply(x, summaryCli,"type",info)
    cgc.table <- do.call(rbind, cgc.table)
    
    #cgc.table <- rbind(cgc.age, cgc.table)
    cgc.table$V3 <- round(as.numeric(cgc.table$V3),digits = 3)
    colnames(cgc.table)[1] <- "characteristics"
    colnames(cgc.table)[4] <- "pvalue"
    #openxlsx::write.xlsx(cgc.table, file = paste0(outdir, 'characteristics.xlsx'), row.names = F)
    return(cgc.table)
    
  }
  
}


