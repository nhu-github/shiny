

# survival
ggsurvival <- function(mx, time = 'OS_months', event = 'OS_status',
                       strata = 'RRAP', showLegend = TRUE, 
                       plottitle = NULL, ytitle = 'Probability of survival', showTable = T){
  require(survival)
  require(survminer)
  #mx <- cgc.mx; time = 'OS_months'; event = 'OS_status'; strata = 'RRAP'
  mx[['time']] <- mx[[time]]
  mx[['event']] <- mx[[event]]
  mx[['strata']] <- mx[[strata]]
  
  if(is.numeric(showLegend)){
    legendPos <- showLegend
  }else if(showLegend){
    legendPos <- c(0.8, 0.9)
  }else{
    legendPos <- 'none'
  }
  
  # legendlab <- c('RRAP-Wt', 'RRAP_Mut')
  # if(is.factor(mx[['strata']])) legendlab <- levels(mx[['strata']])
  
  # if(is.numeric(timelimit)){
  #   mx[mx[['time']] > timelimit, 'event'] <- 0
  #   mx[mx[['time']] > timelimit, 'time'] <- timelimit
  # }
  
  
  mx.cox <- coxph(Surv(time, event)~strata, data = mx)
  mx.cox <- summary(mx.cox)
  mx.cox.HR <- mx.cox$conf.int[c(1,3,4)]
  mx.cox.HR <- round(mx.cox.HR, digits = 2)
  HRlabel <- paste0('HR = ', mx.cox.HR[1], '(', mx.cox.HR[2], ' - ', mx.cox.HR[3], ')')
  mx.cox.logrank <- mx.cox$sctest[3]
  mx.cox.logrank <- ifelse(mx.cox.logrank < 0.001, 
                           sprintf('%.2e', mx.cox.logrank),
                           round(mx.cox.logrank, digits = 3))
  plabel <- paste0('p-value = ', mx.cox.logrank)
  
  label.df <- data.frame(npcx = c(0.05, 0.05), npcy = c(0.2, 0.1), label = c(HRlabel, plabel))
  
  mx.km <- survfit(Surv(time, event)~strata, data = mx, type= "kaplan-meier")
  p <- ggsurvplot(mx.km, data = mx, palette = 'npg', conf.int = F,
                  #xlab = 'Time (Months)', ylab = ytitle,
                  #break.x.by=, break.y.by=0.2,
                  conf.int.style='step', legend.title='', title = plottitle,
                  #legend.labs = legendlab, 
                  legend = legendPos,
                  pval = F, pval.method = F,
                  risk.table = showTable, tables.height = 0.25, 
                  risk.table.title = 'Number at risk', risk.table.y.text = F,
                  censor = T, censor.shape = 4, censor.size = 4,
                  ggtheme = theme(axis.title = element_text(face = 'bold', size = 15, colour = 'black'),
                                  plot.title = element_text(face = 'bold', size = 18, hjust = 0.5, colour = 'black'),
                                  axis.text = element_text(face = 'bold', size = 13, colour = 'black'),
                                  legend.text = element_text(face = 'bold', size = 13, colour = 'black'),
                                  legend.key = element_rect(fill='transparent'),
                                  legend.background = element_blank(),
                                  panel.background=element_rect(fill='transparent'),
                                  axis.line = element_line(size = 0.5)))
 
return(p)
}
  
  
  
  
# ggforest 

pggforest <- function(mx,Multivar,time = 'OS.time', event = 'OS',Stage=NULL){
  require(survival)
  require(forestmodel)
  
  CC <- unique(c(time,event,Stage,Multivar))
  mx <- mx[,CC]
  mx[['time']] <- as.numeric(mx[[time]])
  mx[['event']] <- mx[[event]]
  
  # if(!is.null(Stage)){
  #   mx[['Stage']] <- gsub("A|B|C|D","",mx[[Stage]])
  # }else{
  #   mx[['Stage']] <-  mx[[Stage]]
  # }
  
  panels <- list(list(width = 0.015),
                 list(width = 0.1, display = ~variable, 
                      fontface = "bold", heading = "Variable"),
                 list(width = 0.1, display = ~level),
                 list(width = 0.05, display = ~n, hjust = 1, heading = "N"),
                 list(width = 0.02, item = "vline", hjust = 0.5),
                 list(width = 0.67, item = "forest", hjust = 0.5, heading = "Hazard ratio", 
                      linetype = "dashed", line_x = 0),
                 list(width = 0.1, display = ~ifelse(reference, "Reference", 
                                                     sprintf("%0.2f (%0.2f, %0.2f)",
                                                             trans(estimate), 
                                                             trans(conf.low), 
                                                             trans(conf.high))), 
                      display_na = NA),
                 list(width = 0.02, item = "vline", hjust = 0.5),
                 list(width = 0.04,
                      display = ~ifelse(reference, "", format.pval(p.value, digits = 1, eps = 0.001)),
                      display_na = NA, hjust = 1, heading = "p-value"),
                 list(width = 0.015)
  )
  
  y = Surv(mx[['time']],mx[['event']])
  
  MultivarN <- colnames(mx)[!colnames(mx) %in% c("time","event",time,event)]
  
  fml=as.formula(paste0("y~",paste0(MultivarN,collapse = "+")))
  coxmx=coxph(fml,data = mx)
  p <- forest_model(coxmx, panels, recalculate_width = T, recalculate_height = T,
                    #format_options = list(text_size = 12*textsizescale)
  )
  p$layers[[2]]$aes_params$colour <- 'blue'
  p$layers[[3]]$aes_params$colour <- 'black'
  p$layers[[3]]$aes_params$size <- 0.8
  p <- p + theme(axis.text.x = element_text(size = 12, colour = 'black'))
  return(p)
}

  
  