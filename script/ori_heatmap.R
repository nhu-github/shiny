ori_heatmap <- function(expr,genelist=NULL,cli=NULL,resam=F,show_column_names=T,
                        show_row_names=F,cluster_columns=T,cluster_rows=F){
  require(ComplexHeatmap)
  
  if(!is.null(genelist)){
    expr1 <- as.matrix(expr[rownames(expr)%in%genelist[,1],])
    expr2 <- t(scale(t(expr1)))    
  }else{
    expr1 <- as.matrix(expr[1:50,])
    expr2 <- t(scale(t(expr1))) 
  }
  
  if(!is.null(cli)){
    if(resam==F){
      set.seed(538)
    }else{
      N <- sample(1000,1)
      set.seed(N)
    }
    samp <- intersect(colnames(expr),row.names(cli))
    cli <- cli[rownames(cli)%in%samp,]
    annoInfo <- as.data.frame(apply(as.data.frame(cli), 2, as.factor))
    topAnno <- ComplexHeatmap::HeatmapAnnotation(
      df = annoInfo,
      annotation_legend_param = list(title_gp=gpar(fontsize=8, fontface='bold'), 
                                     labels_gp=gpar(fontsize=8, fontface='bold')),
      annotation_name_gp = gpar(fontsize = 8, fontface='bold'),
      show_legend = T
    )
    
    expr2 <- expr2[,colnames(expr2) %in% samp]
    # expr 
    p <- Heatmap(expr2,name = "RNA",
                 cluster_columns=cluster_columns,
                 cluster_rows=cluster_rows,
                 show_heatmap_legend=T,
                 # heatmap_legend_param = list(legend_height = unit(3, "cm"), 
                 #                           title_position = "topcenter"),
                 show_row_dend=F,
                 show_row_names=show_row_names,
                 show_column_names=show_column_names,
                 top_annotation=topAnno
    )
    
  }else{
    p <- Heatmap(expr2,name = "RNA",
                 cluster_columns=cluster_columns,
                 cluster_rows=cluster_rows,
                 show_heatmap_legend=T,
                 # heatmap_legend_param = list(legend_height = unit(3, "cm"), 
                 #                           title_position = "topcenter"),
                 show_row_dend=F,
                 show_row_names=show_row_names,
                 show_column_names=show_column_names,
                 # top_annotation=topAnno
    )    
    
  }
  
  return(p)
  
}
