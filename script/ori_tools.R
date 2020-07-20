#' @import ggplot2
NULL
#' Filter variants within genes with multiple variants in one patients

#' @param data A dataframe with column ID, gene, and variants, ...
#' @param varorder Priority of variants

rmvar <- function(data, varorder){
  rmg <- function(d, varorder){
    ugene <- unique(d[[2]])
    ud <- lapply(ugene,
                 function(x){
                   udd <- d[d[[2]] == x,];
                   udd <- unique(udd);
                   if(nrow(udd) == 1){
                     return(udd)
                   }else{
                     for(i in varorder){
                       if(i %in% udd[[3]]) return(udd[udd[[3]] == i,])
                     }
                   }
                 })
    ud <- do.call(rbind, ud)
    ud <- as.data.frame(ud, stringsAsFactors = F)
    return(ud)
  }

  uid <- unique(data[[1]])
  udata <- lapply(uid, function(x){rmg(data[data[[1]] == x,], varorder)})
  udata <- do.call(rbind, udata)
  udata <- as.data.frame(udata, stringsAsFactors = F)
  return(udata)
}

## barplot core
# ori_ggbar <- function(ggdata, x, y, group, color){
#   miny <- min(ggdata[[y]], na.rm = T)
#   maxy <- max(ggdata[[y]], na.rm = T)
#   if(miny >= 0 & maxy < 1){
#     ybreak <- pretty(0:(maxy*100))/100
#   }else if(miny >= 0 & maxy > 1){
#     ybreak <- pretty(0:maxy)
#   }else if(miny < 0 & miny >= -1){
#     ybreak <- pretty((miny*100):(maxy*100))/100
#   }else if(miny < -1 & maxy > 1){
#     ybreak <- pretty(miny:maxy)
#   }
#   legendrow <- ceiling(length(levels(ggdata[[group]]))/3)
# 
#   p <- ggplot(ggdata, aes_string(x, y))
#   p <- p + geom_bar(aes_string(fill = group), stat = 'identity',
#                     position = position_stack(), width = 0.5)
#   p <- p + scale_fill_manual(values = color)
#   p <- p + scale_y_continuous(breaks = ybreak, labels = abs(ybreak))
#   p <- p + guides(fill = guide_legend(nrow = legendrow, byrow = T))
#   if(maxy > 1) p <- p + labs(y = 'Counts')
#   if(maxy <= 1) p <- p + labs(y = 'Percentage (%)')
#   p <- p + theme_bw()
#   p <- p + theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_text(size = 8, face = 'bold', color = 'black'),
#     axis.text.y = element_text(size = 8, face = 'bold', color = 'black'),
#     axis.text.x = element_text(size = 8, face = 'bold', color = 'black',
#                                angle = 60, hjust = 1,
#                                margin = margin(2, 0, -4, 0)),
#                                #hjust = 1, vjust = 0.5),
#     plot.title = element_blank(),
#     panel.grid = element_blank(),
#     panel.border = element_blank(),
#     axis.ticks = element_blank(),
#     legend.position = 'bottom',
#     legend.title = element_text(size = 7, face = 'bold', color = 'black'),
#     legend.text = element_text(size = 7, face = 'bold', color = 'black'),
#     legend.key.size = unit(0.5, 'line')
#   )
#   if(miny < 0){
#     p <- p + geom_hline(yintercept = 0, color = 'white', size = 0.3)
#     p <- p + theme(
#       axis.line.y = element_line(color='black', size = 0.3),
#       axis.ticks.y = element_line(color='black', size = 0.3))
#   }else{
#     p <- p + theme(
#       axis.line = element_line(color='black', size = 0.3),
#       axis.ticks = element_line(color='black', size = 0.3))
#   }
# 
#   return(p)
# }

## barplot core v2

ori_ggbar <- function(ggdata, x, y, group, color,ytype){
  miny <- min(ggdata[[y]], na.rm = T)
  #maxy <- max(ggdata[[y]], na.rm = T)
  maxy <- max(unlist(lapply(split(ggdata,ggdata[[x]]), function(x){sum(as.numeric(x[,3]))})),na.rm = T)
  if(miny >= 0 & maxy <= 1){
    ybreak <- pretty(0:(maxy*100))/100
  }else if(miny >= 0 & maxy >= 1){
    ybreak <- pretty(0:maxy)
  }else if(miny < 0 & miny >= -1){
    ybreak <- pretty((miny*100):(maxy*100))/100
  }else if(miny < -1 & maxy >= 1){
    ybreak <- pretty(miny:maxy)
  }
  #legendrow <- ceiling(length(levels(ggdata[[group]]))/3)
  legendrow <- 1
  p <- ggplot(ggdata, aes_string(x, y))
  p <- p + geom_bar(aes_string(fill = group), stat = 'identity',
                    position = position_stack(), width = 0.5)
  p <- p + scale_fill_manual(values = color)
  p <- p + scale_y_continuous(breaks = ybreak, labels = abs(ybreak))
  p <- p + guides(fill = guide_legend(nrow = legendrow, byrow = T))
  
  #if(maxy > 1) p <- p + labs(y = 'Counts')
  #if(maxy <= 1) p <- p + labs(y = 'Percentage (%)')
  if(ytype =='counts') p <- p + labs(y = 'Counts')
  if(ytype == 'percentage') p <- p + labs(y = 'Percentage (%)')
  
  p <- p + theme_bw()
  p <- p + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8, face = 'bold', color = 'black'),
    axis.text.y = element_text(size = 8, face = 'bold', color = 'black'),
    axis.text.x = element_text(size = 8, face = 'bold', color = 'black',
                               angle = 60, hjust = 1,
                               margin = margin(2, 0, -4, 0)),
    #hjust = 1, vjust = 0.5),
    plot.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    legend.position = 'bottom',
    legend.title = element_text(size = 7, face = 'bold', color = 'black'),
    legend.text = element_text(size = 7, face = 'bold', color = 'black'),
    legend.key.size = unit(0.5, 'line')
  )
  if(miny < 0){
    p <- p + geom_hline(yintercept = 0, color = 'white', size = 0.3)
    p <- p + theme(
      axis.line.y = element_line(color='black', size = 0.3),
      axis.ticks.y = element_line(color='black', size = 0.3))
  }else{
    p <- p + theme(
      axis.line = element_line(color='black', size = 0.3),
      axis.ticks = element_line(color='black', size = 0.3))
  }
  return(p)
}







































## profiling
ori_coreoncoprint <- function(mx,
                              rowdf = NULL,
                              topAnno = NULL,
                              color,
                              rowsplit = NULL,
                              colsplit = NULL,
                              cell_width = 0.6,
                              cell_heigh = 0.8){

  # fill cells
  alter_fun = function(x,y,w,h,v){
    n=sum(v)
    h = h*cell_heigh
    w = w*cell_width
    grid.rect(x, y, w, h, gp = gpar(fill = "grey90", col = NA))
    if(n) grid.rect(x,y-h*0.5+1:n/n*h, w, 1/n*h,
                    gp=gpar(fill=color[names(which(v))], col=NA),
                    just = "top")
  }
  # onco
  onco <- oncoPrint(mx,
                    get_type = function(x) unique(strsplit(x, ":")[[1]]),
                    alter_fun = alter_fun,
                    col = color,
                    heatmap_legend_param = list(title='Alternations',
                                                title_position=c('topleft'),
                                                ncol = 1,
                                                title_gp=gpar(fontsize=9, fontface='bold'),
                                                labels_gp=gpar(fontsize=9, fontface='bold')),
                    #row and column title
                    row_title_side='left',
                    row_title_gp=gpar(fontsize=9,
                                      fontface='bold'),
                    column_title_side='bottom',
                    column_title_gp=gpar(col='white'),
                    #row order and column order
                    row_order = rownames(mx),
                    column_order = colnames(mx),
                    #row names
                    show_row_names = T,
                    row_names_side = 'left',
                    row_names_gp = gpar(fontsize = 9,
                                        fontface='bold'),
                    #column names
                    show_column_names = F,
                    column_names_side = 'left',
                    column_names_gp = gpar(fontsize = 9,
                                        fontface='bold'),
                    #


                    top_annotation = topAnno,
                    column_order = colnames(alteration),
                    column_names_gp = gpar(fontsize = 9, fontface='bold'),

                    gap = unit(0.034, 'npc'),
                    show_column_names = F,
                    show_row_barplot = F,
                    show_pct = T,
                    pct_gp=gpar(fontsize=9, col='black', fontface='bold'))

  if(!(is.null(rowdf))){
    rowAnno <- rowAnnotation(bar= row_anno_barplot(rowpercent,
                                                   baseline = 0,
                                                   bar_width = 0.7,
                                                   gp = gpar(fill = color, col = color),
                                                   axis_gp = gpar(fontsize=10, fontface='bold'),
                                                   axis = T, axis_side = 'top', border = F),
                             width = unit(0.09, 'npc'),
                             show_annotation_name = F)
  }

  onco <- onco + rowAnno

  draw(onco, ht_gap = unit(0.01, 'npc'),
       heatmap_legend_side='right', annotation_legend_list = list(Ltype))

}


cal_freq = function(mut, total){
  # 计算基因突变频率
  mut = unique(mut[, c("ORDER_ID", "GENE")])
  mut = as.data.frame(table(mut$GENE))
  colnames(mut) = c("GENE", "MUTANT")
  mut$N = total
  mut$FREQUENCY = mut$MUTANT / mut$N
  mut$WT = mut$N - mut$MUTANT
  mut = mut[order(mut$FREQUENCY, decreasing = T), ]
  return(mut)
}

convert_mut = function(mut=mut, freq, all_sn, n = 20){
  # 将数据转化成画profiling需要的格式
  mut = unique(mut[, c("ORDER_ID", "GENE", "VAR_TYPE_SX")])
  mut = data.table::dcast(mut, ORDER_ID~GENE, fun.aggregate = function(x){paste(x, collapse = ";")})
  #mut <- na.omit(mut)
  rownames(mut) = mut$ORDER_ID
  mut[setdiff(all_sn, mut$ORDER_ID), ] = ""
  mut =data.frame(mut)
  mut = base::t(mut[,-1])
  freq = head(freq, n = n)
  rownames(mut)<- gsub("\\.","-",rownames(mut))
  mut = mut[rownames(mut) %in% freq$GENE, ]
  return(mut)
}

# common mutation of the two groups 
convert_mut1 = function(mut=mut, freq1,freq2,all_sn, n = 20){
  # 将数据转化成画profiling需要的格式
  mut = unique(mut[, c("ORDER_ID", "GENE", "VAR_TYPE_SX")])
  mut = data.table::dcast(mut, ORDER_ID~GENE, fun.aggregate = function(x){paste(x, collapse = ";")})
  #mut <- na.omit(mut)
  rownames(mut) = mut$ORDER_ID
  mut[setdiff(all_sn, mut$ORDER_ID), ] = ""
  mut =data.frame(mut)
  mut = base::t(mut[,-1])
  commongene <- intersect(freq1$GENE,freq2$GENE)
  freq <- freq1[freq1$GENE %in% commongene,]
  freq = head(freq, n = n)
  mut = mut[rownames(mut) %in% freq$GENE, ]
  return(mut)
}


# v1
# ori_HeatmapAnnotation <- function(cli,bar,feature,resam=F,cutoff=NULL){
#   if(resam==F){
#     set.seed(123)
#   }else{
#     N <- sample(1000,1)
#     set.seed(N)
#   }
#   # the cutoff for bar,if NULL,use the median of bar
#   require("grid")
#   annoInfo <- data.frame(cli[,c("ORDER_ID",bar,feature)])
#   rownames(annoInfo) <- annoInfo[,"ORDER_ID"]
#   annoInfo <- data.frame(annoInfo[,-1])
#   annoInfo[[bar]] <- as.numeric(annoInfo[[bar]])
#   annoInfo <- annoInfo[order(annoInfo[[bar]], decreasing = T), ]
#   if(is.null(cutoff)){
#     annoInfo$new <- ifelse(annoInfo[[bar]] > median(annoInfo[[bar]]),"High","Low") 
#   }else{
#     annoInfo$new <- ifelse(annoInfo[[bar]] > cutoff,"High","Low") 
#   }
#   annoInfo <- annoInfo[,c(1,ncol(annoInfo),2:(ncol(annoInfo)-1))]
#   colnames(annoInfo)[2] <- paste0(bar,"1")
#   topAnno <- ComplexHeatmap::HeatmapAnnotation(tmp =anno_barplot(annoInfo[,1], gp = gpar(fill = "#1F78B4"),
#                                                                  height = unit(2, "cm"),width=unit(1, "cm")
#                                                                  ),
#                                                df = annoInfo[,-1:-2],
#                                                annotation_legend_param = list(title_gp=gpar(fontsize=8, fontface='bold'), 
#                                                                               labels_gp=gpar(fontsize=8, fontface='bold')),
#                                                annotation_name_gp = gpar(fontsize = 8, fontface='bold'),
#                                                show_legend = F
#   )
#   topAnno@anno_list$tmp@name <- bar
#   if(length(feature)==1){
#    topAnno@anno_list$df@name <- feature
#   }
#   names(topAnno@anno_list) <- c(bar,feature)
#   
#   return(list("annoInfo"=annoInfo,"topAnno"=topAnno))
# }

# v2
ori_HeatmapAnnotation <- function(annoInfo,feature,bar=NULL,resam=F,annotationsize=8){
  if(resam==F){
    set.seed(123)
  }else{
    N <- sample(1000,1)
    set.seed(N)
  }
  # the cutoff for bar,if NULL,use the median of bar
  require("grid")
  if(!is.null(bar)){
    topAnno <- ComplexHeatmap::HeatmapAnnotation(tmp =anno_barplot(annoInfo[,1], gp = gpar(fill = "#1F78B4"),
                                                                   height = unit(2, "cm"),width=unit(1, "cm")
    ),
    df = annoInfo[,-1:-2],
    annotation_legend_param = list(title_gp=gpar(fontsize=annotationsize, fontface='bold'), 
                                   labels_gp=gpar(fontsize=annotationsize, fontface='bold')),
    annotation_name_gp = gpar(fontsize = annotationsize, fontface='bold'),
    show_legend = F
    )
    topAnno@anno_list$tmp@name <- bar
    if(length(feature)==1){
      topAnno@anno_list$df@name <- feature
    }
    names(topAnno@anno_list) <- c(bar,feature)
    return(list("topAnno"=topAnno))
  }else{
    topAnno <- ComplexHeatmap::HeatmapAnnotation(
      df = annoInfo,
      annotation_legend_param = list(title_gp=gpar(fontsize=annotationsize, fontface='bold'), 
                                     labels_gp=gpar(fontsize=annotationsize, fontface='bold')),
      annotation_name_gp = gpar(fontsize = annotationsize, fontface='bold'),
      show_legend = F
    )
    names(topAnno@anno_list) <- c(feature)
    return(list("topAnno"=topAnno))
  }
  
}



change_mut<- function(mut){
  mutN <- mut
  mutN[mutN==""] = as.numeric(1)
  mutN[!mutN==1] = as.numeric(0)
  mutN <- data.frame(apply(mutN,2,  as.numeric))
  rownames(mutN) <- row.names(mut)
  colnames(mutN) <- colnames(mut)
  mutN <- data.frame(t(mutN))
  return(mutN)
}


ori_HeatmapAnnotation_preorder <- function(cli,bar=NULL,feature,cutoff=NULL){
  require("grid")
  if(!is.null(bar)){
    annoInfo <- data.frame(cli[,c("ORDER_ID",bar,feature)])
    rownames(annoInfo) <- annoInfo[,"ORDER_ID"]
    annoInfo <- data.frame(annoInfo[,-1])
    annoInfo[[bar]] <- as.numeric(annoInfo[[bar]])
    annoInfo <- annoInfo[order(annoInfo[[bar]], decreasing = T), ]
    if(is.null(cutoff)){
      annoInfo$new <- ifelse(annoInfo[[bar]] > median(annoInfo[[bar]]),"High","Low") 
    }else{
      annoInfo$new <- ifelse(annoInfo[[bar]] > cutoff,"High","Low") 
    }
    annoInfo <- annoInfo[,c(1,ncol(annoInfo),2:(ncol(annoInfo)-1))]
    colnames(annoInfo)[2] <- paste0(bar,"1")
    return(annoInfo)
    
  }else{
    annoInfo <- data.frame(cli[,c("ORDER_ID",feature)])
    rownames(annoInfo) <- annoInfo[,"ORDER_ID"]
    annoInfo <- annoInfo[,-1,drop=F]
    return(annoInfo)
  }
  
}




