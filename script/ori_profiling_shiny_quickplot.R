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

plot_landscape = function(mut,cli=NULL,bar,feature,prefix,quickplot=F,waterfall=F,resam=F,cutoff=NULL,n = 30,rownamessize=8){
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
  
  
  
  # alter_fun = function(x,y,w,h,v){
  #   n=sum(v)
  #   h = h*0.8
  #   w = w*0.7
  #   grid.rect(x, y, w, h, gp = gpar(fill = "grey95", col = NA))
  #   if(n) grid.rect(x,y-h*0.5+1:n/n*h, w, 1/n*h, gp=gpar(fill=oncocol[names(which(v))], col=NA), just = "top")
  # }
  
  # alter_fun = list(
  #   background = function(x, y, w, h){
  #     h = h*0.8
  #     w = w*0.7
  #     grid.rect(x, y, w, h, gp = gpar(fill = "grey95", col = NA))
  #   } ,
  #   `Substitution/Indel` = function(x, y, w, h){
  #     h = h*0.8
  #     w = w*0.7
  #     grid.rect(x, y, w, h, gp=gpar(fill="#228B22", col=NA), just = "top")
  #   } ,
  #   `Gene Amplification` = function(x, y, w, h){
  #     h = h*0.8
  #     w = w*0.7
  #     grid.rect(x, y, w, h, gp=gpar(fill="#EE0000", col=NA), just = "top")
  #   } ,
  #   `Gene Homozygous Deletion` = function(x, y, w, h) {
  #     h = h*0.8
  #     w = w*0.7
  #     grid.rect(x, y, w, h, gp=gpar(fill="#0000EE", col=NA), just = "top")
  #   },
  #   
  #   `Fusion/Rearrangement` = function(x, y, w, h){
  #     h = h*0.8
  #     w = w*0.7
  #     grid.rect(x, y, w, h, gp=gpar(fill="#EEEE00", col=NA), just = "top")
  #   },
  #   `Truncation`=function(x, y, w, h){
  #     h = h*0.8
  #     w = w*0.7
  #     grid.rect(x, y, w, h, gp=gpar(fill="#8E388E", col=NA), just = "top")
  #   }
  # )

  varty <- unique(mut[["VAR_TYPE_SX"]])
  if(length(varty)==1){
    alter_fun = function(x,y,w,h,v){
      n=sum(v)
      h = h*0.8
      w = w*0.7
      grid.rect(x, y, w, h, gp = gpar(fill = "grey95", col = NA))
      grid.rect(x,y-h*0.5+1:n/n*h, w, 1/n*h, gp=gpar(fill=oncocol[varty], col=NA), just = "top")
    }
    }else{
      alter_fun = function(x,y,w,h,v){
        n=sum(v)
        h = h*0.8
        w = w*0.7
        grid.rect(x, y, w, h, gp = gpar(fill = "grey95", col = NA))
        if(n) grid.rect(x,y-h*0.5+1:n/n*h, w, 1/n*h, gp=gpar(fill=oncocol[names(which(v))], col=NA), just = "top")
      }
    
  }
  
  varorder = c('Fusion/Rearrangement', 'Substitution/Indel', 
               'Gene Amplification', 'Gene Homozygous Deletion', 
               'Truncation')
  mut <- na.omit(mut[,c("ORDER_ID","GENE","VAR_TYPE_SX")])
  mut <- rmvar(mut,varorder)
  
  allsamples = unique(mut[["ORDER_ID"]])
  
  if(!is.null(cli)){
    HA = ori_HeatmapAnnotation(cli,bar,feature,resam=resam,cutoff=NULL,annotationsize=rownamessize)
    anno_legend_list = lapply(HA$topAnno@anno_list[c(feature)], 
                              function(anno) color_mapping_legend(anno@color_mapping, plot = FALSE))
    for (f in feature) {
      anno_legend_list[[f]]@grob$children[[1]]$gp$fontsize <- rownamessize + 2
      anno_legend_list[[f]]@grob$children[[2]]$children[[1]]$gp$fontsize <- rownamessize + 2
    }

    if(length(feature)==1){
      anno_legend_list[[feature]]@grob$children[[1]][[1]] <- feature 
    }
    
    mut = as.data.frame(mut, stringsAsFactors = F)
    freq = cal_freq(mut, length(allsamples))
    mut = convert_mut(mut, freq, allsamples, n = n)
    mutN <- change_mut(mut)
    mutN <- merge(HA$annoInfo,mutN,by="row.names")
    mutN <- mutN[do.call(order,c(mutN[,-1],list(decreasing=TRUE))),]
    colnames(mutN)[1] <- "ORDER_ID"
    
    mut = mut[, match(mutN$ORDER_ID, colnames(mut))]
    
    #cairo_pdf(paste(prefix, "_landscape.pdf", sep=""), width = 10, height = 6, family = "GB1")
    #pdf(paste(prefix, "_landscape.pdf", sep=""), width = 10, height = 6, onefile = FALSE)
    if(waterfall==F){
      
      if(quickplot==F){
        onco = oncoPrint(mut, get_type = function(x) unique(strsplit(x, ";")[[1]]),
                         alter_fun = alter_fun, col = oncocol,
                         heatmap_legend_param = list(title='Alternations', title_position=c('topleft'), nrow=length(names(table(mut)))-1,
                                                     title_gp=gpar(fontsize=rownamessize+2, fontface='bold'), labels_gp=gpar(fontsize=rownamessize+2)),
                         column_names_gp = gpar(fontsize = rownamessize, fontface='bold'),
                         row_names_gp = gpar(fontsize = rownamessize, fontface='bold'),
                         row_title_gp=gpar(fontsize=rownamessize+2, fontface='bold'),
                         show_column_names = F,
                         show_heatmap_legend=T,
                         show_pct = T,
                         remove_empty_columns = T,
                         top_annotation =HA$topAnno,
                         column_order = mutN$ORDER_ID
        )
        
        # draw(onco, heatmap_legend_side='right',heatmap_legend_list=c(list(lgd_mut),anno_legend_list))
      }else{
        onco <- Heatmap(mut,
                        col = oncocol,
                        na_col = 'grey99',
                        column_names_gp = gpar(fontsize = rownamessize),
                        row_names_gp = gpar(fontsize = rownamessize, fontface='bold'),
                        #row_names_side='left',
                        #row_title_side='left',
                        heatmap_legend_param = list(title='Alternations', title_position=c('topleft'), nrow=length(names(table(mut)))-1,
                                                    title_gp=gpar(fontsize=rownamessize+2, fontface='bold'), labels_gp=gpar(fontsize=rownamessize+2)),
                        row_title_gp=gpar(fontsize=rownamessize,fontface='bold'),
                        top_annotation = HA$topAnno,
                        column_order = mutN$ORDER_ID,
                        #row_order = colnames(mutN1),
                        row_order = as.character(freq$GENE[1:n]),
                        show_column_names = F,
                        show_heatmap_legend = T,
                        use_raster = T,
                        raster_device = "png",
                        raster_quality = 2)
        #draw(onco, heatmap_legend_side='right',heatmap_legend_list=c(list(lgd_mut),anno_legend_list))
        
      }
      
      
    }else{
      if(quickplot==F){
        onco = oncoPrint(mut, get_type = function(x) unique(strsplit(x, ";")[[1]]),
                         alter_fun = alter_fun, col = oncocol,
                         heatmap_legend_param = list(title='Alternations', title_position=c('topleft'), nrow=length(names(table(mut)))-1,
                                                     title_gp=gpar(fontsize=rownamessize+2, fontface='bold'), labels_gp=gpar(fontsize=rownamessize+2)),
                         column_names_gp = gpar(fontsize = rownamessize, fontface='bold'),
                         row_names_gp = gpar(fontsize = rownamessize, fontface='bold'),
                         row_title_gp=gpar(fontsize=rownamessize, fontface='bold'),
                         show_column_names = F,
                         show_heatmap_legend=T,
                         show_pct = T,
                         remove_empty_columns = T,
                         top_annotation =HA$topAnno
                         #  column_order = rownames(HA$annoInfo)
        )
        # draw(onco, heatmap_legend_side='right',heatmap_legend_list=c(list(lgd_mut),anno_legend_list))
        
      }else{
        
        mutNN <- change_mut(mut)
        mutNN <- mutNN[,as.character(freq$GENE[1:n])]
        mutNN <- mutNN[do.call(order,mutNN),]
        # colnames(mutN)[1] <- "ORDER_ID"
        
        onco <- Heatmap(mut,
                        col = oncocol,
                        na_col = 'grey99',
                        column_names_gp = gpar(fontsize = rownamessize),
                        row_names_gp = gpar(fontsize = rownamessize, fontface='bold'),
                        #row_names_side='left',
                        #row_title_side='left',
                        heatmap_legend_param = list(title='Alternations', title_position=c('topleft'), nrow=length(names(table(mut)))-1,
                                                    title_gp=gpar(fontsize=rownamessize+2, fontface='bold'), labels_gp=gpar(fontsize=rownamessize+2)),
                        row_title_gp=gpar(fontsize=rownamessize,fontface='bold'),
                        top_annotation = HA$topAnno,
                        column_order =  rownames(mutNN),
                        #row_order = colnames(mutN1),
                        row_order = as.character(freq$GENE[1:n]),
                        show_column_names = F,
                        show_heatmap_legend = T,
                        use_raster = T,
                        raster_device = "png",
                        raster_quality = 2)
        #draw(onco, heatmap_legend_side='right',heatmap_legend_list=c(list(lgd_mut),anno_legend_list))
        
      }
      
    }
    
    p = draw(onco, heatmap_legend_side='right',heatmap_legend_list=c(anno_legend_list))
    print(p)
  }else{
    mut = as.data.frame(mut, stringsAsFactors = F)
    freq = cal_freq(mut, length(allsamples))
    mut = convert_mut(mut, freq, allsamples, n = n)
    
    if(quickplot==F){
      onco = oncoPrint(mut, get_type = function(x) unique(strsplit(x, ";")[[1]]),
                       alter_fun = alter_fun, col = oncocol,
                       heatmap_legend_param = list(title='Alternations', title_position=c('topleft'), nrow=length(names(table(mut)))-1,
                                                   title_gp=gpar(fontsize=rownamessize+2, fontface='bold'), labels_gp=gpar(fontsize=rownamessize+2)),
                       column_names_gp = gpar(fontsize = rownamessize, fontface='bold'),
                       row_names_gp = gpar(fontsize = rownamessize, fontface='bold'),
                       row_title_gp=gpar(fontsize=rownamessize, fontface='bold'),
                       show_column_names = F,
                       show_heatmap_legend=T,
                       show_pct = T,
                       remove_empty_columns = T)
      # draw(onco, heatmap_legend_side='right',heatmap_legend_list=c(list(lgd_mut)))
      
    }else{
      
      mutNN <- change_mut(mut)
      mutNN <- mutNN[,as.character(freq$GENE[1:n])]
      mutNN <- mutNN[do.call(order,mutNN),]
      # colnames(mutN)[1] <- "ORDER_ID"
      
      onco <- Heatmap(mut,
                      col = oncocol,
                      na_col = 'grey99',
                      column_names_gp = gpar(fontsize =rownamessize),
                      row_names_gp = gpar(fontsize = rownamessize, fontface='bold'),
                      #row_names_side='left',
                      #row_title_side='left',
                      heatmap_legend_param = list(title='Alternations', title_position=c('topleft'), nrow=length(names(table(mut)))-1,
                                                  title_gp=gpar(fontsize=rownamessize+2, fontface='bold'), labels_gp=gpar(fontsize=rownamessize+2)),
                      row_title_gp=gpar(fontsize=rownamessize,fontface='bold'),
                      # top_annotation = HA$topAnno,
                      column_order =  rownames(mutNN),
                      #row_order = colnames(mutN1),
                      row_order = as.character(freq$GENE[1:n]),
                      show_column_names = F,
                      show_heatmap_legend = T,
                      use_raster = T,
                      raster_device = "png",
                      raster_quality = 2)
      #draw(onco, heatmap_legend_side='right',heatmap_legend_list=c(list(lgd_mut)))
      
    }
    
   # p = draw(onco, heatmap_legend_side='right',heatmap_legend_list=c(list(lgd_mut)))
    p = draw(onco, heatmap_legend_side='right')
    print(p)
    
  }
  
  # dev.off()
  return(p)
}



plot_pathway_cli = function(mut, cli=NULL, bar, feature, prefix, 
                            pathway_genes, quickplot=F , waterfall=F, resam=F, nfreq=0.01, cutoff=NULL,rownamessize=8){
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
  
  varty <- unique(mut[["VAR_TYPE_SX"]])
  if(length(varty)==1){
    alter_fun = function(x,y,w,h,v){
      n=sum(v)
      h = h*0.8
      w = w*0.7
      grid.rect(x, y, w, h, gp = gpar(fill = "grey95", col = NA))
      grid.rect(x,y-h*0.5+1:n/n*h, w, 1/n*h, gp=gpar(fill=oncocol[varty], col=NA), just = "top")
    }
  }else{
    alter_fun = function(x,y,w,h,v){
      n=sum(v)
      h = h*0.8
      w = w*0.7
      grid.rect(x, y, w, h, gp = gpar(fill = "grey95", col = NA))
      if(n) grid.rect(x,y-h*0.5+1:n/n*h, w, 1/n*h, gp=gpar(fill=oncocol[names(which(v))], col=NA), just = "top")
    }
    
  }
  

  varorder = c('Fusion/Rearrangement', 'Substitution/Indel', 
               'Gene Amplification', 'Gene Homozygous Deletion', 
               'Truncation')
  mut <- na.omit(mut[,c("ORDER_ID","GENE","VAR_TYPE_SX")])
  mut <- rmvar(mut,varorder)
  

  if(!is.null(cli)){
    HA = ori_HeatmapAnnotation(cli,bar,feature,resam=resam,cutoff=NULL,annotationsize=rownamessize)
    anno_legend_list = lapply(HA$topAnno@anno_list[c(feature)], 
                              function(anno) color_mapping_legend(anno@color_mapping, plot = FALSE))
    for (f in feature) {
      anno_legend_list[[f]]@grob$children[[1]]$gp$fontsize <- rownamessize + 2
      anno_legend_list[[f]]@grob$children[[2]]$children[[1]]$gp$fontsize <- rownamessize + 2
    }
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
    # add cli,Note the purpose of this step
    mut[setdiff(rownames(HA$annoInfo), mut$ORDER_ID), ] = ""
    
    mut = base::t(mut[,-1])
    mut =  mut[, match(rownames(HA$annoInfo), colnames(mut))]
    
    all_genes = all_genes[match(rownames(mut), all_genes$Gene), ]
    rownames(all_genes) = all_genes$Gene
    all_genes = all_genes[, c("Pathway"), drop = F]
    
    # all_genes_order
    all_genes_order <- merge(all_genes,freq[,c(1,4),drop=F],by.x = "row.names",by.y = "GENE")
    colnames(all_genes_order)[1] <- "GENE"
    all_genes_order <- all_genes_order[do.call(order,c(all_genes_order[,-1],list(decreasing=TRUE))),]
    
    # new column order
    mutN <- change_mut(mut)
    mutN <- merge(HA$annoInfo,mutN,by="row.names")
    mutN <- mutN[do.call(order,c(mutN[,-1],list(decreasing=TRUE))),]
    colnames(mutN)[1] <- "ORDER_ID"
    
    #pdf(paste(prefix, "_pathway_landscape.pdf", sep=""), width = 12, height = 10, onefile = FALSE)
    #cairo_pdf(paste(prefix, "_pathway_landscape.pdf", sep=""), width = 12, height = 10, family = "GB1")
    
    if(waterfall==F){
      if(quickplot==F){
        onco = oncoPrint(mat = mut, alter_fun = alter_fun, 
                         get_type = function(x) unique(strsplit(x, ";")[[1]]),
                         col = oncocol,
                         heatmap_legend_param = list(title='Alternations', title_position=c('topleft'), nrow=length(names(table(mut)))-1,
                                                     title_gp=gpar(fontsize=rownamessize+2, fontface='bold'), labels_gp=gpar(fontsize=rownamessize+2)),
                         column_names_gp =  gpar(fontsize = rownamessize, fontface='bold'),
                         row_names_gp = gpar(fontsize = rownamessize, fontface='bold'),
                         row_title_gp= gpar(fontsize = rownamessize, fontface='bold'),
                         show_pct = T, 
                         row_names_side = "right",
                         show_heatmap_legend=T,
                         gap = unit(3, 'mm'),
                         row_split = all_genes,
                         remove_empty_columns = T,
                         top_annotation = HA$topAnno,
                         column_order = mutN$ORDER_ID,
                         row_order = all_genes_order$GENE
        )
        
        p = draw(onco, heatmap_legend_side='right',heatmap_legend_list=c(anno_legend_list))
        print(p)    
        
      }else{
        onco <- Heatmap(mut,
                        col = oncocol,
                        na_col = 'grey99',
                        column_names_gp = gpar(fontsize = rownamessize, fontface='bold'),
                        row_names_gp = gpar(fontsize = rownamessize, fontface='bold'),
                        #row_names_side='left',
                        #row_title_side='left',
                        heatmap_legend_param = list(title='Alternations', title_position=c('topleft'), nrow=length(names(table(mut)))-1,
                                                    title_gp=gpar(fontsize=rownamessize+2, fontface='bold'), labels_gp=gpar(fontsize=rownamessize+2)),
                        row_title_gp=gpar(fontsize = rownamessize, fontface='bold'),
                        top_annotation = HA$topAnno,
                        column_order = mutN$ORDER_ID,
                        #column_gap = unit(3, "mm"),
                        gap = unit(3, 'mm'),
                        #cell_fun = cell_fun,
                        #row_order = colnames(mutN1),
                        row_order = all_genes_order$GENE,
                        row_split = all_genes,
                        show_column_names = F,
                        show_heatmap_legend = T,
                        use_raster = T,
                        raster_device = "png",
                        raster_quality = 2)
        draw(onco, heatmap_legend_side='right',heatmap_legend_list=c(anno_legend_list))
      }
      
      
    }else{
      
      
      mutNN <- change_mut(mut)
      mutNN <- mutNN[,all_genes_order$GENE]
      mutNN <- mutNN[do.call(order,mutNN),]
      # colnames(mutN)[1] <- "ORDER_ID"
      if(quickplot==F){
        onco = oncoPrint(mat = mut, alter_fun = alter_fun, 
                         get_type = function(x) unique(strsplit(x, ";")[[1]]),
                         col = oncocol,
                         heatmap_legend_param = list(title='Alternations', title_position=c('topleft'), nrow=length(names(table(mut)))-1,
                                                     title_gp=gpar(fontsize=rownamessize+2, fontface='bold'), labels_gp=gpar(fontsize=rownamessize+2)),
                         column_names_gp =gpar(fontsize = rownamessize, fontface='bold'),
                         row_names_gp = gpar(fontsize = rownamessize, fontface='bold'),
                         row_title_gp=gpar(fontsize = rownamessize, fontface='bold'),
                         show_pct = T, 
                         row_names_side = "right",
                         show_heatmap_legend=T,
                         gap = unit(3, 'mm'),
                         row_split = all_genes,
                         remove_empty_columns = T,
                         top_annotation = HA$topAnno,
                         column_order = rownames(mutNN)
        )
        
        p = draw(onco, heatmap_legend_side='right',heatmap_legend_list=c(anno_legend_list))
        print(p)    
        
      }else{
        onco <- Heatmap(mut,
                        col = oncocol,
                        na_col = 'grey99',
                        column_names_gp = gpar(fontsize = rownamessize, fontface='bold'),
                        row_names_gp = gpar(fontsize = rownamessize, fontface='bold'),
                        #row_names_side='left',
                        #row_title_side='left',
                        heatmap_legend_param = list(title='Alternations', title_position=c('topleft'), nrow=length(names(table(mut)))-1,
                                                    title_gp=gpar(fontsize=rownamessize+2, fontface='bold'), labels_gp=gpar(fontsize=rownamessize+2)),
                        row_title_gp=gpar(fontsize = rownamessize, fontface='bold'),
                        top_annotation = HA$topAnno,
                        column_order =  rownames(mutNN),
                        #row_order = colnames(mutN1),
                        row_order = all_genes_order$GENE,
                        show_column_names = F,
                        show_heatmap_legend = T,
                        row_split = all_genes,
                        use_raster = T,
                        raster_device = "png",
                        raster_quality = 2)
        p =draw(onco, heatmap_legend_side='right',heatmap_legend_list=c(anno_legend_list))
        print(p) 
        
      }
      
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
    
    #all_genes_order
    all_genes_order <- merge(all_genes,freq[,c(1,4),drop=F],by.x = "row.names",by.y = "GENE")
    colnames(all_genes_order)[1] <- "GENE"
    all_genes_order <- all_genes_order[do.call(order,c(all_genes_order[,-1],list(decreasing=TRUE))),]
    
    # new order
    mutNN <- change_mut(mut)
    mutNN <- mutNN[,all_genes_order$GENE]
    mutNN <- mutNN[do.call(order,mutNN),]
    
    if(quickplot==F){
      onco = oncoPrint(mat = mut, alter_fun = alter_fun, 
                       get_type = function(x) unique(strsplit(x, ";")[[1]]),
                       col = oncocol,
                       heatmap_legend_param = list(title='Alternations', title_position=c('topleft'), nrow=length(names(table(mut)))-1,
                                                   title_gp=gpar(fontsize=rownamessize+2, fontface='bold'), labels_gp=gpar(fontsize=rownamessize+2)),
                       column_names_gp = gpar(fontsize = rownamessize, fontface='bold'),
                       row_names_gp = gpar(fontsize = rownamessize, fontface='bold'),
                       row_title_gp=gpar(fontsize = rownamessize, fontface='bold'),
                       show_pct = T, 
                       row_names_side = "right",
                       show_heatmap_legend=T,
                       gap = unit(3, 'mm'),
                       row_split = all_genes,
                       remove_empty_columns = T,
                       #top_annotation = HA$topAnno,
                       row_order = all_genes_order$GENE,
                       column_order = rownames(mutNN)
      )
      #p = draw(onco, heatmap_legend_side='right',heatmap_legend_list=c(list(lgd_mut)))
      p = draw(onco, heatmap_legend_side='right')
      print(p) 
      
    }else{
      onco <- Heatmap(matrix=mut,
                      col = oncocol,
                      na_col = 'grey99',
                      column_names_gp = gpar(fontsize = rownamessize, fontface='bold'),
                      row_names_gp = gpar(fontsize = rownamessize, fontface='bold'),
                      #row_names_side='left',
                      #row_title_side='left',
                      heatmap_legend_param = list(title='Alternations', title_position=c('topleft'), nrow=length(names(table(mut)))-1,
                                                  title_gp=gpar(fontsize=rownamessize+2, fontface='bold'), labels_gp=gpar(fontsize=rownamessize+2)),
                      row_title_gp=gpar(fontsize = rownamessize, fontface='bold'),
                      #top_annotation = HA$topAnno,
                      column_order =  rownames(mutNN),
                      #row_order = colnames(mutN1),
                      row_order = all_genes_order$GENE,
                      show_column_names = F,
                      show_heatmap_legend = T,
                      row_split = all_genes,
                      gap = unit(3, 'mm'),
                      use_raster = T,
                      raster_device = "png",
                      raster_quality = 2)
      
     # p = draw(onco, heatmap_legend_side='right',heatmap_legend_list=c(list(lgd_mut)))
      p = draw(onco, heatmap_legend_side='right')
      print(p) 
      
    }
    
  }
  #dev.off()
  return(p)
}
















