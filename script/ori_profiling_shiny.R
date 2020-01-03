#' @include CountVariants.R ori_tools.R
#' @import ComplexHeatmap
#' @import ggpubr
#' @import Cairo
NULL
#' Drawing two profiling ,one is mutation landscape,one is with pathway 
#'
#' @inheritParams CountVariants
#' @param mut Mutation file, there must be one colnames named "ORDER_ID"
#' @param cli Clinical file
#' @param bar Column names of clinical file used for draw barplot on the top of profiling.Continous variable.
#' @param feature Column names of clinical file used for pheatmap on the top of profiling.Discrete variable.
#' @param prefix output pdf file's prefix
#' @param cutoff The cutoff of groupby.Default is the median of groupby value.
#' @param n  The number of first top mutantdgenes to show.Default is 30.
#' @param all_sn All samples' names
#' @param pathway Pathway dataframe, column as pathway name,rownames as the gene names in the pathway
#' @param nfreq  The cutoff of mutation frequency to show.
#' 
#' @return output pdf file
#' @examples
#' plot_landscape(Mut0,info0,bar = 'TMB',feature=c("GENDER","Smoker"),prefix="test",n = 30)
#' plot_pathway_cli(Mut0, info0,"age","gender","test2", unique(Mut0$ORDER_ID), pathway)
#'
#' @export


plot_landscape = function(mut,cli=NULL,bar,feature,prefix,waterfall=F,resam=F,cutoff=NULL,n = 30){
  require("ComplexHeatmap")
  oncocol <- c("Substitution/Indel" = "#228B22",        ##green
               "Gene Amplification" = "#EE0000",        ##red
               "Gene Homozygous Deletion" = "#0000EE",  ##blue
               "Fusion/Rearrangement" = "#EEEE00",      ##yellow
               "Truncation" = "#8E388E")                ##purple
  #'Splicing' = '#B2DF8A')                 ##
  
  lgd_mut <- Legend(at = c("0","1","2","3","4"),
                    labels = c("Substitution/Indel",
                               "Gene Amplification",
                               "Gene Homozygous Deletion",
                               "Fusion/Rearrangement",
                               "Truncation"),
                    title = "Alternations",
                    legend_gp = gpar(fill = c("#228B22",
                                              "#EE0000",
                                              "#0000EE",
                                              "#EEEE00",
                                              "#8E388E")),
                    labels_gp = gpar( fontsize = 10),
                    title_gp = gpar( fontsize = 10, fontface="bold"))
  
  
  
  alter_fun = function(x,y,w,h,v){
    n=sum(v)
    h = h*0.8
    w = w*0.7
    grid.rect(x, y, w, h, gp = gpar(fill = "grey95", col = NA))
    if(n) grid.rect(x,y-h*0.5+1:n/n*h, w, 1/n*h, gp=gpar(fill=oncocol[names(which(v))], col=NA), just = "top")
  }
  
  varorder = c('Fusion/Rearrangement', 'Substitution/Indel', 
               'Gene Amplification', 'Gene Homozygous Deletion', 
               'Truncation')
  mut <- na.omit(mut[,c("ORDER_ID","GENE","VAR_TYPE_SX")])
  
  mut <- rmvar(mut,varorder)
  
  allsamples = unique(mut[["ORDER_ID"]])
  
  
  if(!is.null(cli)){
    HA = ori_HeatmapAnnotation(cli,bar,feature,resam=resam,cutoff=NULL)
    anno_legend_list = lapply(HA$topAnno@anno_list[c(feature)], 
                              function(anno) color_mapping_legend(anno@color_mapping, plot = FALSE))
    if(length(feature)==1){
      anno_legend_list[[feature]]@grob$children[[1]][[1]] <- feature 
    }
    
    mut = as.data.frame(mut, stringsAsFactors = F)
    freq = cal_freq(mut, length(allsamples))
    mut = convert_mut(mut, freq, allsamples, n = n)
    mut = mut[, match(rownames(HA$annoInfo), colnames(mut))]
    
    #cairo_pdf(paste(prefix, "_landscape.pdf", sep=""), width = 10, height = 6, family = "GB1")
    #pdf(paste(prefix, "_landscape.pdf", sep=""), width = 10, height = 6, onefile = FALSE)
    if(waterfall==F){
      onco = oncoPrint(mut, get_type = function(x) unique(strsplit(x, ";")[[1]]),
                       alter_fun = alter_fun, col = oncocol,
                       heatmap_legend_param = list(title='Alternations', title_position=c('topleft'), nrow=6,
                                                   title_gp=gpar(fontsize=8, fontface='bold'), labels_gp=gpar(fontsize=8, fontface='bold')),
                       column_names_gp = gpar(fontsize = 8, fontface='bold'),
                       row_names_gp = gpar(fontsize = 8, fontface='bold'),
                       row_title_gp=gpar(fontsize=10, fontface='bold'),
                       show_column_names = F,
                       show_heatmap_legend=F,
                       show_pct = T,
                       remove_empty_columns = T,
                       top_annotation =HA$topAnno,
                       column_order = rownames(HA$annoInfo)
      )
    }else{
      onco = oncoPrint(mut, get_type = function(x) unique(strsplit(x, ";")[[1]]),
                       alter_fun = alter_fun, col = oncocol,
                       heatmap_legend_param = list(title='Alternations', title_position=c('topleft'), nrow=6,
                                                   title_gp=gpar(fontsize=8, fontface='bold'), labels_gp=gpar(fontsize=8, fontface='bold')),
                       column_names_gp = gpar(fontsize = 8, fontface='bold'),
                       row_names_gp = gpar(fontsize = 8, fontface='bold'),
                       row_title_gp=gpar(fontsize=10, fontface='bold'),
                       show_column_names = F,
                       show_heatmap_legend=F,
                       show_pct = T,
                       remove_empty_columns = T,
                       top_annotation =HA$topAnno
                       #  column_order = rownames(HA$annoInfo)
      )
      
    }
    
    p = draw(onco, heatmap_legend_side='right',heatmap_legend_list=c(list(lgd_mut),anno_legend_list))
    print(p)
  }else{
    mut = as.data.frame(mut, stringsAsFactors = F)
    freq = cal_freq(mut, length(allsamples))
    mut = convert_mut(mut, freq, allsamples, n = n)
    onco = oncoPrint(mut, get_type = function(x) unique(strsplit(x, ";")[[1]]),
                     alter_fun = alter_fun, col = oncocol,
                     heatmap_legend_param = list(title='Alternations', title_position=c('topleft'), nrow=6,
                                                 title_gp=gpar(fontsize=8, fontface='bold'), labels_gp=gpar(fontsize=8, fontface='bold')),
                     column_names_gp = gpar(fontsize = 8, fontface='bold'),
                     row_names_gp = gpar(fontsize = 8, fontface='bold'),
                     row_title_gp=gpar(fontsize=10, fontface='bold'),
                     show_column_names = F,
                     show_heatmap_legend=F,
                     show_pct = T,
                     remove_empty_columns = T)
    p = draw(onco, heatmap_legend_side='right',heatmap_legend_list=c(list(lgd_mut)))
    print(p)

  }

 # dev.off()
  return(p)
}



plot_pathway_cli = function(mut, cli=NULL, bar, feature, prefix, 
                            pathway_genes,waterfall=F, resam=F, nfreq=0.03, cutoff=NULL){
  ## cli
  oncocol <- c("Substitution/Indel" = "#228B22",        ##green
               "Gene Amplification" = "#EE0000",        ##red
               "Gene Homozygous Deletion" = "#0000EE",  ##blue
               "Fusion/Rearrangement" = "#EEEE00",      ##yellow
               "Truncation" = "#8E388E")                ##purple
  #'Splicing' = '#B2DF8A')                 ##
  
  lgd_mut <- Legend(at = c("0","1","2","3","4"),
                    labels = c("Substitution/Indel",
                               "Gene Amplification",
                               "Gene Homozygous Deletion",
                               "Fusion/Rearrangement",
                               "Truncation"),
                    title = "Alternations",
                    legend_gp = gpar(fill = c("#228B22",
                                              "#EE0000",
                                              "#0000EE",
                                              "#EEEE00",
                                              "#8E388E")),
                    labels_gp = gpar( fontsize = 10),
                    title_gp = gpar( fontsize = 10, fontface="bold"))
  
  
  
  alter_fun = function(x,y,w,h,v){
    n=sum(v)
    h = h*0.8
    w = w*0.8
    grid.rect(x, y, w, h, gp = gpar(fill = "grey95", col = NA))
    if(n) grid.rect(x,y-h*0.5+1:n/n*h, w, 1/n*h, gp=gpar(fill=oncocol[names(which(v))], col=NA), just = "top")
  }
  
  varorder = c('Fusion/Rearrangement', 'Substitution/Indel', 
               'Gene Amplification', 'Gene Homozygous Deletion', 
               'Truncation')
  mut <- na.omit(mut[,c("ORDER_ID","GENE","VAR_TYPE_SX")])
  mut <- rmvar(mut,varorder)
  
  
  
  if(!is.null(cli)){
    HA = ori_HeatmapAnnotation(cli,bar,feature,resam=resam,cutoff=NULL)
    anno_legend_list = lapply(HA$topAnno@anno_list[c(feature)], 
                              function(anno) color_mapping_legend(anno@color_mapping, plot = FALSE))
    if(length(feature)==1){
      anno_legend_list[[feature]]@grob$children[[1]][[1]] <- feature 
    }

    all_sn = unique(mut[["ORDER_ID"]])
    ## pathway 
    pathway_genes = as.data.frame(pathway_genes, stringsAsFactors = F)
    all_genes =  data.table::melt(pathway_genes, measure.vars = colnames(pathway_genes), 
                                  variable.name = "Pathway", value.name = "Gene")
    all_genes = all_genes[!is.na(all_genes$Gene), ]
    freq = cal_freq(mut, length(all_sn))
    highgene <- freq$GENE[freq$FREQUENCY>=nfreq]
    all_genes = all_genes[all_genes$Gene %in% highgene,]
    
    ## Mutation
    mut = as.data.frame(mut, stringsAsFactors = F)
    mut = mut[mut$ORDER_ID %in% all_sn, ]
    mut = unique(mut[mut$GENE %in% all_genes$Gene, c("ORDER_ID", "GENE", "VAR_TYPE_SX")])
    mut = data.table::dcast(mut, ORDER_ID~GENE, fun.aggregate = function(x){paste(x, collapse = ";")})
    rownames(mut) = mut$ORDER_ID
    mut[setdiff(all_sn, mut$ORDER_ID), ] = ""
    mut = base::t(mut[,-1])
    mut =  mut[, match(rownames(HA$annoInfo), colnames(mut))]
    
    all_genes = all_genes[match(rownames(mut), all_genes$Gene), ]
    rownames(all_genes) = all_genes$Gene
    all_genes = all_genes[, c("Pathway"), drop = F]
    
    
    #pdf(paste(prefix, "_pathway_landscape.pdf", sep=""), width = 12, height = 10, onefile = FALSE)
    #cairo_pdf(paste(prefix, "_pathway_landscape.pdf", sep=""), width = 12, height = 10, family = "GB1")
    
    if(waterfall==F){
      onco = oncoPrint(mat = mut, alter_fun = alter_fun, 
                       get_type = function(x) unique(strsplit(x, ";")[[1]]),
                       col = oncocol,
                       heatmap_legend_param = list(title='Alternations', title_position=c('leftcenter'), 
                                                   nrow=1, title_gp=gpar(fontface='bold', size = 10), 
                                                   labels_gp=gpar(fontface='bold', size = 8)),
                       column_names_gp = gpar(fontface='bold'),
                       row_names_gp = gpar(fontface='bold', cex = 0.5),
                       row_title_gp=gpar(fontface='bold'),
                       show_pct = T, 
                       row_names_side = "right",
                       show_heatmap_legend=F,
                       gap = unit(3, 'mm'),
                       row_split = all_genes,
                       remove_empty_columns = T,
                       top_annotation = HA$topAnno,
                       column_order = rownames(HA$annoInfo))
      
      p = draw(onco, heatmap_legend_side='right',heatmap_legend_list=c(list(lgd_mut),anno_legend_list))
      print(p)    
      
    }else{
      onco = oncoPrint(mat = mut, alter_fun = alter_fun, 
                       get_type = function(x) unique(strsplit(x, ";")[[1]]),
                       col = oncocol,
                       heatmap_legend_param = list(title='Alternations', title_position=c('leftcenter'), 
                                                   nrow=1, title_gp=gpar(fontface='bold', size = 10), 
                                                   labels_gp=gpar(fontface='bold', size = 8)),
                       column_names_gp = gpar(fontface='bold'),
                       row_names_gp = gpar(fontface='bold', cex = 0.5),
                       row_title_gp=gpar(fontface='bold'),
                       show_pct = T, 
                       row_names_side = "right",
                       show_heatmap_legend=F,
                       gap = unit(3, 'mm'),
                       row_split = all_genes,
                       remove_empty_columns = T,
                       top_annotation = HA$topAnno,
                       #column_order = rownames(HA$annoInfo)
      )
      
      p = draw(onco, heatmap_legend_side='right',heatmap_legend_list=c(list(lgd_mut),anno_legend_list))
      print(p)    

  }

  }else{
    
    all_sn = unique(mut[["ORDER_ID"]])
    ## pathway 
    pathway_genes = as.data.frame(pathway_genes, stringsAsFactors = F)
    all_genes =  data.table::melt(pathway_genes, measure.vars = colnames(pathway_genes), 
                                  variable.name = "Pathway", value.name = "Gene")
    all_genes = all_genes[!is.na(all_genes$Gene), ]
    freq = cal_freq(mut, length(all_sn))
    highgene <- freq$GENE[freq$FREQUENCY>=nfreq]
    all_genes = all_genes[all_genes$Gene %in% highgene,]
    
    ## Mutation
    mut = as.data.frame(mut, stringsAsFactors = F)
    mut = mut[mut$ORDER_ID %in% all_sn, ]
    mut = unique(mut[mut$GENE %in% all_genes$Gene, c("ORDER_ID", "GENE", "VAR_TYPE_SX")])
    mut = data.table::dcast(mut, ORDER_ID~GENE, fun.aggregate = function(x){paste(x, collapse = ";")})
    rownames(mut) = mut$ORDER_ID
    mut[setdiff(all_sn, mut$ORDER_ID), ] = ""
    mut = base::t(mut[,-1])
    #mut =  mut[, match(rownames(HA$annoInfo), colnames(mut))]
    
    all_genes = all_genes[match(rownames(mut), all_genes$Gene), ]
    rownames(all_genes) = all_genes$Gene
    all_genes = all_genes[, c("Pathway"), drop = F]
    
    onco = oncoPrint(mat = mut, alter_fun = alter_fun, 
                     get_type = function(x) unique(strsplit(x, ";")[[1]]),
                     col = oncocol,
                     heatmap_legend_param = list(title='Alternations', title_position=c('leftcenter'), 
                                                 nrow=1, title_gp=gpar(fontface='bold', size = 10), 
                                                 labels_gp=gpar(fontface='bold', size = 8)),
                     column_names_gp = gpar(fontface='bold'),
                     row_names_gp = gpar(fontface='bold', cex = 0.5),
                     row_title_gp=gpar(fontface='bold'),
                     show_pct = T, 
                     row_names_side = "right",
                     show_heatmap_legend=F,
                     gap = unit(3, 'mm'),
                     row_split = all_genes,
                     remove_empty_columns = T
                     #top_annotation = HA$topAnno,
                     #column_order = rownames(HA$annoInfo)
    )
    p = draw(onco, heatmap_legend_side='right',heatmap_legend_list=c(list(lgd_mut)))
    print(p)    
    
  }
  #dev.off()
  return(p)
}

















