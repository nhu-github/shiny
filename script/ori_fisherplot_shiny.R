#' @param data The result after fisher exact test
#' @param OR  The OR column
#' @param padj  The padj column
#' @param GENE  The column with recording gene names
#' @param outdir output directory
#' 
#' @return A fisher test plot 
#' @examples
#' ori_fisherplot(l,OR="OR",padj ="fisher.test.Padj", GENE="GENE",outdir = "tmp")
#'
#' @export

ori_fisherplot<- function(data,OR,padj,GENE,outdir=NULL){
  require(ggplot2)
  if(!is.null(outdir)){
    outdir <- paste0("./",outdir,"/")
    if(!dir.exists(outdir)){
      dir.create(outdir)
    }
  }
  l <- data 
  l$x <- log2(l[[OR]])
  l[is.infinite(l$x) & l$x > 0,"x"] <- 8
  l[is.infinite(l$x) & l$x < 0,"x"] <- -8
  l$y <- -log2(l[[padj]])
  l$GENE <- l[[GENE]]
  l$padj <- l[[padj]]
  
  p <- ggplot(data = l, mapping = aes(x = x, y = y))+
    geom_point(position = position_jitter())+
    geom_vline(xintercept = 0, linetype = 2)+
    geom_hline(yintercept = min(-log2(0.05), max(l[['y']])), linetype = 2)+
    labs(x = 'log2(OR)', y = '-log2(padj)')+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 12, face = 'bold'),
      axis.title = element_text(size = 14, face = 'bold')
    )
  
  p <- p + ggrepel::geom_text_repel(data = dplyr::filter(l, padj < 0.05), 
                                    aes(label = GENE), nudge_y = 0.2)
  #ggsave(paste0(outdir,'GENEvs', '_fisher.pdf'), height = 8, width = 8)
  return(p)
}

# l <- t$fisher_test
# l$GENE <- row.names(l)
# p <- ori_fisherplot(l,OR="OR",padj ="fisher.test.Padj", GENE="GENE",outdir = "tmp")







