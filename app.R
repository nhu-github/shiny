#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
options(shiny.maxRequestSize=70*1024^2)
library(shiny)
library(openxlsx)
library(ggplot2)
library(shinythemes)
library(DT)

#source("./script/ori_profiling_shiny.R")
#source("./script/ori_profiling_shiny_quickplot_origin.R")
#source("./script/ori_profiling_shiny_quickplot.R")
source("./script/ori_profiling_shiny_quickplot1.R")
source("./script/ori_boxplot_shiny.R")
source("./script/ori_barplot.R")
source("./script/CountVariants.R")
source("./script/ori_tools.R")
source("./script/lollipop.R")
source("./script/ori_lollipop_shiny.R")
source("./script/ori_statcli_shiny.R")
source("./script/ori_fisherplot_shiny.R")
#source("./script/ori_biomarker_gene_shiny.R")
source("./script/ori_biomarker_gene_shiny1.R")
source("./script/ori_survival.R")
source("./script/gene_ex_co_shiny.r")
source("./script/drive_gene_shiny.r")
source("./script/circus_shiny.r")
source("./script/CNV_shiny.r")
source("./script/ori_mechine_learning.R")
source("./script/xCell.R")
source("./script/ori_scattercor_plot.R")
source("./script/ori_heatmap.R")
source("./script/CIBERSORT.R")
source("./script/ori_tcgaCompare.R")

# Define UI for application 
ui <- fluidPage(
  
  theme = shinytheme('cerulean'),
    #shinythemes::themeSelector()
  #theme = shinytheme('sandstone'),
    navbarPage("Himalaya 1.1", 
               
               tabPanel("Mutation",
                        fluidRow(
                            column(2,
                                   wellPanel(
                                       h4(strong("Input your mutation file")),
                                       selectInput("file1",label= "choose an example or your own data", 
                                                   choices = c("Example"="Example1", "Your own data" = "load_my_own1")),
                                       conditionalPanel("input.file1 == 'Example1'",
                                                        downloadButton('downloadEx1', 'Download example')),
                                       conditionalPanel("input.file1 == 'load_my_own1'",
                                                        fileInput('file31', 'Choose xlsx File', 
                                                                  accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv', '.txt')))
                                       
                                   ),
                            
                                   wellPanel(
                                       h4(strong("Input your clinical file")),
                                       selectInput("file2",label= "choose an example or your own data", 
                                                   choices = c("Example"="Example2", "Your own data" = "load_my_own2")),
                                       conditionalPanel("input.file2 == 'Example2'",
                                                        downloadButton('downloadEx2', 'Download example')),
                                       conditionalPanel("input.file2 == 'load_my_own2'",
                                                        fileInput('file32', 'Choose xlsx File', 
                                                                  accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv', '.txt')))
                                       
                                   ),
                                   
                                   conditionalPanel("input.cPanels1 == 2",
                                                    wellPanel(
                                                      h4(strong("Profiling")),
                                                      
                                                      selectInput("select01", "Divided into two subgroups", 
                                                                  choices=c("TMB"),multiple = F),
                                                      selectInput("select011", "Continuous or Discrete ",
                                                                  choices=c("continuous","discrete"),
                                                                  selected = "continuous",
                                                                  multiple = FALSE),
                                                      selectInput("select11", "Choose clinical feature", 
                                                                  choices=c("GENDER"),multiple = T),
                                                      sliderInput("number", 
                                                                  label = "The number of first top mutant genes:",
                                                                  min = 10, max = 100, value = 30, step = 1
                                                      ),
                                                      sliderInput("rownamessize",
                                                                  label = "The size of genenames",
                                                                  min = 6, max = 16, value = 8, step = 1
                                                      ),
                                                      
                                                      selectInput("selectwaterfall", "show the waterfall", 
                                                                  choices=c("YES"=TRUE,"NO"=FALSE),multiple = F,
                                                                  selected = FALSE),
                                                      
                                                      selectInput("selectquickplot", "quick plot", 
                                                                  choices=c("YES"=TRUE,"NO"=FALSE),multiple = F,
                                                                  selected = FALSE),
                                                      
                                                      selectInput("selectchangecolor", "change the color of legend", 
                                                                  choices=c("YES"=TRUE,"NO"=FALSE),multiple = F,
                                                                  selected = FALSE),
                                                      actionButton("refresh", "Refresh"),
                                                      selectInput("select21", "Show the pathway", 
                                                                  choices = c("NO" = "NO","YES"="YES")
                                                      ),
                                                      conditionalPanel("input.select21 == 'YES'",
                                                                       selectInput("file3",label= "choose an example or your own data", 
                                                                                   choices = c("Example"="Example3", "Your own data" = "load_my_own3")),
                                                                       conditionalPanel("input.file3 == 'Example3'",
                                                                                        downloadButton('downloadEx3', 'Download example')),
                                                                       conditionalPanel("input.file3 == 'load_my_own3'",
                                                                                        fileInput('file33', 'Choose xlsx File', 
                                                                                                  accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv', '.txt'))
                                                                                        
                                                                       )
                                                      ),
                                                      
                                                      #radioButtons("pathway", "Show the pathway", choices = c("Yes" = "Yes", "No" = "NO"), selected = "NO"),
                                                      
                                                      
                                                      sliderInput("nfreq", 
                                                                  label = "The mutation frequency of signal pathway",
                                                                  min = 0, max = 0.1, value = 0.02, step = 0.01
                                                      )
                                                      #textInput("prefix", "filename"),
                                                      
                                                      # hr()                                
                                                    )),
                                   
                                   conditionalPanel("input.cPanels1 == 2",
                                                    h4(strong("Download")),
                                                    wellPanel(
                                                        textInput("fname22", "filename", value = "profiling"),
                                                        downloadButton('Downloadprofiling', 'Download profiling')
                                                    )),
                                   
                                   conditionalPanel("input.cPanels1 == 3",
                                                    wellPanel(
                                                        h4(strong("Boxplot")),
                                                        selectInput("selectbox", "select feature", 
                                                                    choices=c("TMB")),
                                                        selectInput("select_gene", "Gene", 
                                                                    c("TP53"), choices=c("TP53"), multiple = F),
                                                        selectInput("select_type", "choose", 
                                                                    choices=c("boxplot","violin"),
                                                                    selected = "boxplot",
                                                                    multiple = F),
                                                        hr()                                
                                                    )),
                                   
                                   conditionalPanel("input.cPanels1 == 3",
                                                    h4(strong("Download")),
                                                    wellPanel(
                                                        textInput("fname23", "filename", value = "boxplot"),
                                                        downloadButton('Downloadbox', 'Download boxplot')
                                                    )),
                                   
                                   conditionalPanel("input.cPanels1 == 4",
                                                    wellPanel(
                                                        h4(strong("Barplot")),
                                                        selectInput("selectbarplotyaxis", "show the y axis", 
                                                                    choices = c("counts","percentage"),
                                                                    selected = "percentage"
                                                        ),
                                                        selectInput("selectbarplot", "Divided into two subgroup", 
                                                                    choices = c("NO" = "NO","YES"="YES"),
                                                                    selected = "NO"
                                                        ),
                                                        conditionalPanel("input.selectbarplot == 'YES'",
                                                                         selectInput("selectBarfeature",label= "select feature", 
                                                                                     choices = c("TMB")),
                                                                         selectInput("selectBargroupty",label= "Continuous or Discrete", 
                                                                                     choices = c("continuous","discrete"),
                                                                                     selected = "continuous",
                                                                                     multiple = FALSE),
                                                                         selectInput("selectBarbyorder",label= "Reorder the barplot", 
                                                                                     choices=c("YES"=TRUE,"NO"=FALSE),
                                                                                     multiple = F,
                                                                                     selected = FALSE)
                                                        ),
                                                        sliderInput("numberbarplotcutoff", 
                                                                    label = "The value of cutoff:",
                                                                    min = 5, max = 50, value = 10, step = 1
                                                        ),
 
                                                        sliderInput("numberbarplot", 
                                                                    label = "The number of first top mutant genes:",
                                                                    min = 10, max = 100, value = 30, step = 1
                                                        ),
                                                        hr()                                
                                                    )),
                             
                                   conditionalPanel("input.cPanels1 == 4",
                                                    h4(strong("Download")),
                                                    wellPanel(
                                                        textInput("fname24", "filename", value = "barplot"),
                                                        downloadButton('Downloadbar', 'Download barplot')
                                                    )),
                                                     
                                   conditionalPanel("input.cPanels1 == 5",
                                                    wellPanel(
                                                        h4(strong("Choose genes")),
                                                        selectInput("select51", "choose one gene", 
                                                                    c("TP53"), choices=c("TP53"), multiple = FALSE),
                                                        selectInput("select52", "refseq ID", 
                                                                    choices=c("refSeqID"), multiple = FALSE),
                                                        hr()                                
                                                    )),
                                   
                                   conditionalPanel("input.cPanels1 == 5",
                                                    h4(strong("Download")),
                                                    wellPanel(
                                                        textInput("fname25", "filename", value = "lillipop"),
                                                        downloadButton('Downloadlollipop', 'Download loillipop')
                                                    )),
                                   
                                   conditionalPanel("input.cPanels1 == 6",
                                                    wellPanel(
                                                      h4(strong("Three mutation types")),
                                                      selectInput("select61", "choose", 
                                                                  choices=c("Gene mutation frequency"="GeneFre",
                                                                            "Variation frequency" = "VarFre",
                                                                            "Pathway mutation frequency" = "PathwayFre"
                                                                  ),
                                                                  multiple = FALSE),
                                                      
                                                      selectInput("selectstat", "Divided into two subgroup", 
                                                                  choices = c("NO" = "NO","YES"="YES"),
                                                                  selected = "NO"
                                                      ),
                                                      
                                                      conditionalPanel("input.selectstat == 'YES'",
                                                                       selectInput("selectstatfeaturetype",label= "Continuous or Discrete",
                                                                                   choices = c("Continuous","Discrete"),
                                                                                   selected = "Continuous"
                                                                       ),
                                                                       selectInput("selectstatfeature",label= "select feature", 
                                                                                   choices = c("TMB"))
                                                      ),
                                                      
                                                      
                                                      sliderInput("numbercountvariantcutoff", 
                                                                  label = "The value of cutoff:",
                                                                  min = 5, max = 50, value = 12, step = 1
                                                      ),
                                                      
                                                    )),
                                   
                                   conditionalPanel("input.cPanels1 == 6",
                                                    h4(strong("Download")),
                                                    wellPanel(
                                                      textInput("fname26", "filename", value = "CountVariant"),
                                                      downloadButton('Downloadstat', 'Download stat')
                                                    )),
                                   
                                   conditionalPanel("input.cPanels1 == 7",
                                                    wellPanel(
                                                      h4(strong("Clinical statistics")),
                                                      selectInput("select71", "Divided into two subgroups",
                                                                  choices=c("TMB"),
                                                                  multiple = FALSE),
                                                      selectInput("select72", "Continuous or Discrete ",
                                                                  choices=c("continuous","discrete"),
                                                                  selected = "continuous",
                                                                  multiple = FALSE),
                                                      sliderInput("numberclicutoff", 
                                                                  label = "The value of cutoff:",
                                                                  min = 5, max = 50, value = 12, step = 1
                                                      ),
                                                      selectInput("selectclifea", "Select clinical features",
                                                                  choices = c("Age"),
                                                                  selected = "Age",
                                                                  multiple = TRUE
                                                      )

                                                    )),
                               

                                   conditionalPanel("input.cPanels1 == 7",
                                                    h4(strong("Download")),
                                                    wellPanel(
                                                      textInput("fname27", "filename", value = "ClinicalStatistics"),
                                                      downloadButton('Downloadclistat', 'Download clinical')
                                                    )),
                                   
                                   conditionalPanel("input.cPanels1 == 8",
                                                    wellPanel(
                                                      h4(strong("Statistical test")),
                                                      selectInput("select80", "Continuous or Discrete ",
                                                                  choices=c("continuous","discrete"),
                                                                  selected = "continuous",
                                                                  multiple = FALSE),
                                                      
                                                      selectInput("select81", "Divided into two subgroups",
                                                                  choices=c("TMB"),
                                                                  multiple = FALSE),
                                                      
                                                      sliderInput("numbertestcutoff", 
                                                                  label = "The value of cutoff:",
                                                                  min = 5, max = 100, value = 12, step = 1
                                                      ),
                                                      
                                                      sliderInput("selectmutfreq", 
                                                                  label = "Gene mutation frequency:",
                                                                  min = 0, max = 100, value = 5, step = 1
                                                      ),
                                                      selectInput("select82", "choose statistical test type", 
                                                                  choices=c("Wilcoxon test and t test"="Wiltest",
                                                                            "fisher test and chisq test" = "fishertest"
                                                                  ),
                                                                  multiple = FALSE),

                                                    )),
                                   
                                   conditionalPanel("input.cPanels1 == 8",
                                                    h4(strong("Download")),
                                                    wellPanel(
                                                      textInput("fname28", "filename", value = "StatisticalTest"),
                                                      downloadButton('Downloadtest', 'Download Statistical Test')
                                                    )),
                                   
                                   conditionalPanel("input.cPanels1 == 9",
                                                    wellPanel(
                                                      h4(strong("fisher test plot")),
                                                      selectInput("select90", "Continuous or Discrete ",
                                                                  choices=c("continuous","discrete"),
                                                                  selected = "continuous",
                                                                  multiple = FALSE),
                                                      selectInput("select91", "Divided into two subgroups",
                                                                  choices=c("TMB"),
                                                                  multiple = FALSE),
                                                      sliderInput("numberfisherteplotcutoff", 
                                                                  label = "The value of cutoff:",
                                                                  min = 5, max = 100, value = 12, step = 1
                                                      ),
                                                      sliderInput("selectmutfreq9", 
                                                                  label = "Gene mutation frequency:",
                                                                  min = 0, max = 100, value = 5, step = 1
                                                      ),
                                                    )),
                                   
                                   conditionalPanel("input.cPanels1 == 9",
                                                    h4(strong("Download")),
                                                    wellPanel(
                                                      textInput("fname29", "filename", value = "fishertest_plot"),
                                                      downloadButton('Downloadfisherplot', 'Download fishertest plot')
                                                    )),
                                   
                                   #10 
                                   conditionalPanel("input.cPanels1 == 10",
                                                    wellPanel(
                                                      h4(strong("Co-mutation plot")),
                                                      selectInput("select101", "gene colnames",
                                                                  choices=c("GENE"),
                                                                  multiple = FALSE),
                                                      
                                                      selectInput("select102", "sample colnames",
                                                                  choices=c("ORDER_ID"),
                                                                  multiple = FALSE),
                                                      
                                                      selectInput("select103", "VAR TYPE colnames",
                                                                  choices=c("VAR_TYPE"),
                                                                  multiple = FALSE),
                                                      
                                                      sliderInput("select104", 
                                                                  label = "gene names size",
                                                                  min = 0.1, max = 1.5, value = 0.8, step = 0.1
                                                      ),

                                                      sliderInput("numberCoMutation", 
                                                                  label = "Number of genes",
                                                                  min = 5, max = 50, value = 30, step = 1
                                                      ),

                                                    )),
                                   
                                   conditionalPanel("input.cPanels1 == 10",
                                                    h4(strong("Download")),
                                                    wellPanel(
                                                      textInput("fname210", "filename", value = "Co-mutation plot"),
                                                      downloadButton('DownloadCoMutationplot', 'Download Co-mutation plot')
                                                    )),
                                   
                                   #find driver gene
                                   conditionalPanel("input.cPanels1 == 11",
                                                    wellPanel(
                                                      h4(strong("Driver gene plot")),
                                                      selectInput("select111", "gene colnames",
                                                                  choices=c("GENE"),
                                                                  multiple = FALSE),
                                                      selectInput("select112", "sample colnames",
                                                                  choices=c("ORDER_ID"),
                                                                  multiple = FALSE),
                                                      selectInput("select113", "VAR TYPE colnames",
                                                                  choices=c("VAR_TYPE"),
                                                                  multiple = FALSE),
                                                      selectInput("select114", "AA CHANGE colnames",
                                                                  choices=c("AA_CHANGE"),
                                                                  multiple = FALSE),

                                                    )),
                                   
                                   conditionalPanel("input.cPanels1 == 11",
                                                    h4(strong("Download")),
                                                    wellPanel(
                                                      textInput("fname211", "filename", value = "Driver gene"),
                                                      downloadButton('DownloadDriverGeneplot', 'Download driver gene plot'),
                                                      br(),
                                                      br(),
                                                      downloadButton('DownloadDriverGenetable', 'Download driver gene table')
                                                      )),
                                   
                                   # circos
                                   conditionalPanel("input.cPanels1 == 12",
                                                    wellPanel(
                                                      h4(strong("Driver gene plot")),
                                                      selectInput("selectMutType", "Select mutation type",
                                                                  choices=c("snv","cnv","fusion"),
                                                                  selected =c("snv","cnv","fusion"),
                                                                  multiple = T),
                                           
                                                      selectInput("select121", "gene colnames",
                                                                  choices=c("GENE"),
                                                                  multiple = FALSE),
                                                      selectInput("select122", "var type colnames",
                                                                  choices=c("VAR_TYPE"),
                                                                  multiple = FALSE),
                                                      selectInput("select123", "genomic colnames",
                                                                  choices=c("GENOMIC"),
                                                                  multiple = FALSE),
                                                      selectInput("select124", "gene pair colnames",
                                                                  choices=c("DNA_CHANGE"),
                                                                  multiple = FALSE),
                                                      selectInput("select125", "amp del colnames",
                                                                  choices=c("VAR_TYPE_SX"),
                                                                  multiple = FALSE),
                                                      
                                                    )),
                                   
                                   conditionalPanel("input.cPanels1 == 12",
                                                    h4(strong("Download")),
                                                    wellPanel(
                                                      textInput("fname212", "filename", value = "Circos"),
                                                      downloadButton('DownloadCircosplot', 'Download Circos plot'),
                                                    )),
                                   
                                   # CNV
                                   conditionalPanel("input.cPanels1 == 13",
                                                    wellPanel(
                                                      h4(strong("Driver gene plot")),
                                                      selectInput("selectVariantType", "Select variant type",
                                                                  choices=c("CNV"),
                                                                  selected =c("CNV"),
                                                                  multiple = F),
                                                      selectInput("select131", "gene colnames",
                                                                  choices=c("GENE"),
                                                                  multiple = FALSE),
                                                      selectInput("select132", "sample colnames",
                                                                  choices=c("ORDER_ID"),
                                                                  multiple = FALSE),
                                                      selectInput("select133", "cnv type colnames",
                                                                  choices=c("DNA_CHANGE"),
                                                                  multiple = FALSE),
                                                      selectInput("select134", "genomic colnames",
                                                                  choices=c("GENOMIC"),
                                                                  multiple = FALSE),
                                                      sliderInput("select135", 
                                                                  label = "The number of genes with the highest frequency ",
                                                                  min = 1, max = 30, value = 10, step = 1
                                                      ),
                                                      sliderInput("select136", 
                                                                  label = "The size of gene name",
                                                                  min = 0, max = 2, value = 0.6, step = 0.1
                                                      ),
                                                      sliderInput("select137", 
                                                                  label = "The size of cytoband ",
                                                                  min = 0, max = 2, value = 0.6, step = 0.1
                                                      ),
                                                      
                                                    )),
                                   
                                   conditionalPanel("input.cPanels1 == 13",
                                                    h4(strong("Download")),
                                                    wellPanel(
                                                      textInput("fname213", "filename", value = "CNV"),
                                                      downloadButton('DownloadCNVplot', 'Download CNV plot'),
                                                      br(),
                                                      br(),
                                                      downloadButton('DownloadCNVtable', 'Download CNV table')
                                                    )),
                                                      
                                   conditionalPanel("input.cPanels1 == 14",
                                                    wellPanel(
                                                      h4(strong("TMBcompare")),
                                                      textInput("selectcohortName", "cohortName",
                                                                value=c("ori") ),
                                                      selectInput("selectmedianCol", "medianCol",
                                                                  choices=c("red","gray70"),
                                                                  selected=c("red"),
                                                                  multiple = F),
                                                      selectInput("selectTCGAcohort", "TCGA cohort color",
                                                                  choices=c("gray70","black"),
                                                                  selected = c("gray70"),
                                                                  multiple = F),
                                                      
                                                      selectInput("selectInputcohort", "Input cohort color",
                                                                  choices=c("gray70","black"),
                                                                  selected = c("black"),
                                                                  multiple = F),
                                                      
                                                      

                                                    )),
                                   
                                   conditionalPanel("input.cPanels1 == 14",
                                                    h4(strong("Download")),
                                                    wellPanel(
                                                      textInput("fnameTMBcompare", "filename", value = "TMBcompare"),
                                                      downloadButton('DownloadTMBcompareplot', 'Download TMBcompare plot'),
                                                      br(),
                                                      br(),
                                                      downloadButton('DownloadTMBcomparetable', 'Download TMBcompare table')
                                                    )),
                                   

                            ),
                                   
                                  
                            column(10,
                                   tabsetPanel(
                                       tabPanel("Manual", htmlOutput("ReadMe1"), value =1),
                                       tabPanel("Profiling", htmlOutput("pv22"), plotOutput("profiling", height= 1000, width = 1500), value = 2),
                                       tabPanel("Boxplot", htmlOutput("pv31"), plotOutput("boxplot", height= 800, width = 800), htmlOutput("pv32"), value =3),
                                       tabPanel("Barplot", htmlOutput("pv41"), plotOutput("barplot", height= 800, width = 1000), value =4),
                                       tabPanel("Lollipop", htmlOutput("pv25"), plotOutput("lollipop", height= 600,width = 1000), value = 5),
                                       tabPanel("Count variant", htmlOutput("pv61"), DT::dataTableOutput("stat",width = 1200), value = 6),
                                       tabPanel("Clinical statistics", htmlOutput("pv71"), tableOutput("clinicalStat"), value = 7),
                                       tabPanel("Statistical test", htmlOutput("pv81"), DT::dataTableOutput("statisticaltest",width = 1200), value = 8),
                                       tabPanel("Fishertest plot", htmlOutput("pv91"), plotOutput("fishertestplot",height= 800, width = 1000), value = 9),
                                       tabPanel("Co-mutation plot", htmlOutput("pv101"), plotOutput("CoMutationplot",height= 800, width = 1000), value = 10),
                                       tabPanel("Find driver gene", htmlOutput("pv111"), plotOutput("drivergeneplot", height= 800, width = 1000), htmlOutput("pv112"), DT::dataTableOutput("drivergenetable",width = 800), value =11),
                                       tabPanel("Circos", htmlOutput("pv121"), plotOutput("circosplot",height= 1200, width = 1200), value = 12),
                                       tabPanel("CNV", htmlOutput("pv131"), plotOutput("CNVplot",height= 800, width = 1000), htmlOutput("pv132"), DT::dataTableOutput("CNVtable",width = 800), value =13),
                                       tabPanel("TMBcompare", htmlOutput("pvsession222"), plotOutput("TMBcompareplot", height= 800, width = 1000),DT::dataTableOutput("TMBcomparetable",width = 800), value = 14),
                                       id = "cPanels1"
                                   )                
                                   
                            ),
                            
                            column(12,
                                   tags$head(tags$style(type="text/css", "
                                                    #loadmessage {
                                                    position: fixed;
                                                    bottom: 0px;
                                                    right: 0px;
                                                    width: 100%;
                                                    padding: 5px 0px 5px 0px;
                                                    text-align: center; 
                                                    font-weight: bold;
                                                    font-size: 100%;
                                                    color: #000000;
                                                    background-color: #b8b8b8;
                                                    z-index: 105;
                                                    }
                                                    ")),
                                   conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                    tags$div("Loading...",id="loadmessage"))
                                   
                            )
                            
                          
                            
                        )
               ),
               
               tabPanel("Survival",
                        fluidRow(
                          column(2,
                                 # wellPanel(
                                 #   h4(strong("Input your mutation file")),
                                 #   selectInput("file1",label= "choose an example or your own data",
                                 #               choices = c("Example"="Example1", "Your own data" = "load_my_own1")),
                                 #   conditionalPanel("input.file1 == 'Example1'",
                                 #                    downloadButton('downloadEx1', 'Download example')),
                                 #   conditionalPanel("input.file1 == 'load_my_own1'",
                                 #                    fileInput('file31', 'Choose xlsx File',
                                 #                              accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv', '.txt')))
                                 #
                                 # ),
                                 #
                                 wellPanel(
                                   h4(strong("Input your survival file")),
                                   selectInput("filesur2",label= "choose an example or your own data",
                                               choices = c("Example"="Examplesur2", "Your own data" = "load_my_ownsur2")),
                                   conditionalPanel("input.filesur2 == 'Examplesur2'",
                                                    downloadButton('downloadExsur2', 'Download example')),
                                   conditionalPanel("input.filesur2 == 'load_my_ownsur2'",
                                                    fileInput('filesur32', 'Choose xlsx File',
                                                              accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv', '.txt')))

                                 ),

                                 conditionalPanel("input.cPanelssession1 == 2",
                                                  wellPanel(
                                                    h4(strong("Survival")),
                                                    selectInput("selecttime", "Time",
                                                                choices=c("OS.time"),multiple = F),
                                                    selectInput("selectevent", "Event",
                                                                choices=c("OS"),multiple = F),

                                                    selectInput("selectstrata", "Strata",
                                                                choices=c("Mutation"),multiple = F),

                                                    selectInput("selectshowTable", "Show the table",
                                                                choices=c("YES"=TRUE,"NO"=FALSE),multiple = F,
                                                                selected = TRUE),
                                                    


                                                  )),

                                 conditionalPanel("input.cPanelssession1 == 2",
                                                  h4(strong("Download")),
                                                  wellPanel(
                                                    textInput("fnamesurvival", "filename", value = "survival"),
                                                    downloadButton('Downloadsurvival', 'Download survival')
                                                  )),

                                 conditionalPanel("input.cPanelssession1 == 3",
                                                  wellPanel(
                                                    h4(strong("ggforest")),
                                                    selectInput("selectforesttime", "Time",
                                                                choices=c("OS.time"),multiple = F),
                                                    selectInput("selectggevent", "Event",
                                                                choices=c("OS"),multiple = F),
                                                    # selectInput("selectstage", "Stage",
                                                    #             choices=c("pathologic_stage"),multiple = F),
                                                    selectInput("selectMultivariate", "Multivariate",
                                                                choices=c("Age", "gender"),
                                                                selected = c("Age", "gender"),
                                                                multiple = T),
                                                    hr()
                                                  )),

                                 conditionalPanel("input.cPanelssession1 == 3",
                                                  h4(strong("Download")),
                                                  wellPanel(
                                                    textInput("fnamegggforest", "filename", value = "ggforest"),
                                                    downloadButton('Downloadggforest', 'Download ggforest')
                                                  )),





                          ),


                          column(10,
                                 tabsetPanel(
                                   tabPanel("Manual", htmlOutput("ReadMe2"), value =1),
                                   tabPanel("Survival", htmlOutput("pvsession22"), plotOutput("survivalplot", height= 800, width = 1000), value = 2),
                                   tabPanel("ggforest", htmlOutput("pvsession31"), plotOutput("ggforest", height= 800, width = 800), htmlOutput("pvsession32"), value =3),
                                   # tabPanel("Barplot", htmlOutput("pv41"), plotOutput("barplot", height= 800, width = 1000), value =4),
                                   # tabPanel("Lollipop", htmlOutput("pv25"), plotOutput("lollipop", height= 600,width = 1000), value = 5),
                                   # tabPanel("Count variant", htmlOutput("pv61"), DT::dataTableOutput("stat",width = 1200), value = 6),
                                   # tabPanel("Clinical statistics", htmlOutput("pv71"), tableOutput("clinicalStat"), value = 7),
                                   # tabPanel("Statistical test", htmlOutput("pv81"), DT::dataTableOutput("statisticaltest",width = 1200), value = 8),
                                   # tabPanel("Fishertest plot", htmlOutput("pv91"), plotOutput("fishertestplot",height= 800, width = 1000), value = 9),

                                   id = "cPanelssession1"
                                 )

                          ),

                          # column(12,
                          #        tags$head(tags$style(type="text/css", "
                          #                           #loadmessage {
                          #                           position: fixed;
                          #                           bottom: 0px;
                          #                           right: 0px;
                          #                           width: 100%;
                          #                           padding: 5px 0px 5px 0px;
                          #                           text-align: center;
                          #                           font-weight: bold;
                          #                           font-size: 100%;
                          #                           color: #000000;
                          #                           background-color: #b8b8b8;
                          #                           z-index: 105;
                          #                           }
                          #                           ")),
                          #        conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                          #                         tags$div("Loading...",id="loadmessage"))
                          #
                          # )



                        )
               ),
               
               
               
               tabPanel("Machine Learning",
                        fluidRow(
                          column(2,
                                 # wellPanel(
                                 #   h4(strong("Input your mutation file")),
                                 #   selectInput("file1",label= "choose an example or your own data",
                                 #               choices = c("Example"="Example1", "Your own data" = "load_my_own1")),
                                 #   conditionalPanel("input.file1 == 'Example1'",
                                 #                    downloadButton('downloadEx1', 'Download example')),
                                 #   conditionalPanel("input.file1 == 'load_my_own1'",
                                 #                    fileInput('file31', 'Choose xlsx File',
                                 #                              accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv', '.txt')))
                                 #
                                 # ),
                                 #
                                 wellPanel(
                                   h4(strong("Input your file")),
                                   selectInput("fileml2",label= "choose an example or your own data",
                                               choices = c("Example"="Exampleml2", "Your own data" = "load_my_ownml2")),
                                   conditionalPanel("input.fileml2 == 'Exampleml2'",
                                                    downloadButton('downloadExml2', 'Download example')),
                                   conditionalPanel("input.fileml2 == 'load_my_ownml2'",
                                                    fileInput('fileml32', 'Choose xlsx File',
                                                              accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv', '.txt')))
                                   
                                 ),
                                 
                                 conditionalPanel("input.cPanelssession2 == 2",
                                                  wellPanel(
                                                    h4(strong("Machine learning")),
                                                    selectInput("selectmethod", "Method",
                                                                choices=c("Random Forest"="rf",
                                                                          "Support Vector Machine"="svmLinear",
                                                                          "Logistic Regression"="glm",
                                                                          "Naive Bayes"="nb",
                                                                          "k-Nearest Neighbors"="knn",
                                                                          "Decision Tree"="rpart"
                                                                ),
                                                                multiple = F),
                                                    selectInput("selectLabelcol", "label colname",
                                                                choices=c("diabetes"),multiple = F),
                                                    
                                                    sliderInput("selectenCV", 
                                                                label = "number of cross validation",
                                                                min = 2, max = 20, value = 5, step = 1),
                                                    
                                                    sliderInput("selecteSplitData", 
                                                                label = "Split data into training and testing sets",
                                                                min = 0, max = 1, value = 0.8, step = 0.1),
                                                    
                                                    
                                                    
                                                  )),
                                 
                                 
                                 conditionalPanel("input.cPanelssession2 == 2",
                                                  h4(strong("Download")),
                                                  wellPanel(
                                                    textInput("fnameml", "filename", value = "MachineLearning"),
                                                    downloadButton('DownloadMLplot', 'Download ML plot'),
                                                    br(),
                                                    br(),
                                                    downloadButton('DownloadMLimportplot', 'Download ML Importantce plot')
                                                  )),
                                 
                                 
                          ),
                          
                          
                          column(10,
                                 tabsetPanel(
                                   tabPanel("Manual", htmlOutput("ReadMe3"), value =1),
                                   tabPanel("MechineLearning", htmlOutput("pvsessionML"), plotOutput("MLcmplot", height= 600, width =800), plotOutput("MLimportplot", height= 600, width = 800),value = 2),
                                   
                                   id = "cPanelssession2"
                                 )
                          ),
                        )
               ),
               tabPanel("Expression",
                        fluidRow(
                          column(2,
                                 # wellPanel(
                                 #   h4(strong("Input your mutation file")),
                                 #   selectInput("file1",label= "choose an example or your own data",
                                 #               choices = c("Example"="Example1", "Your own data" = "load_my_own1")),
                                 #   conditionalPanel("input.file1 == 'Example1'",
                                 #                    downloadButton('downloadEx1', 'Download example')),
                                 #   conditionalPanel("input.file1 == 'load_my_own1'",
                                 #                    fileInput('file31', 'Choose xlsx File',
                                 #                              accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv', '.txt')))
                                 #
                                 # ),
                                 #
                                 wellPanel(
                                   h4(strong("Input your file")),
                                   selectInput("fileex2",label= "choose an example or your own data",
                                               choices = c("Example"="Exampleex2", "Your own data" = "load_my_ownex2")),
                                   conditionalPanel("input.fileex2 == 'Exampleex2'",
                                                    downloadButton('downloadExex2', 'Download example')),
                                   conditionalPanel("input.fileex2 == 'load_my_ownex2'",
                                                    fileInput('fileex32', 'Choose xlsx File',
                                                              accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv', '.txt')))
                                   
                                 ),
                                 
                                 #xCell
                                 conditionalPanel("input.cPanelssession3 == 2",
                                                  wellPanel(
                                                    h4(strong("xCell")),
                                                    h5("This may take some time"),
                                                    selectInput("selectsig", "Choose gene signatures",
                                                                choices=c("xCell(N=64)"),multiple = F),
                                                    actionButton("startxCell","start")
                                                  )),


                                 conditionalPanel("input.cPanelssession3 == 2",
                                                  h4(strong("Download")),
                                                  wellPanel(
                                                    textInput("fnamexcell", "filename", value = "xCell"),
                                                    downloadButton('Downloadxcell', 'Download xCell table'),

                                                  )),
                                 
                                 # cibersort
                                 conditionalPanel("input.cPanelssession3 == 3",
                                                  wellPanel(
                                                    h4(strong("Cibersort")),
                                                    h5("This may take some time"),
                                                    selectInput("selectmixture_file", "Choose mixture file",
                                                                choices=c("LM22"),multiple = F),
                                                    actionButton("startcibersort","start")
                                                  )),

                                 conditionalPanel("input.cPanelssession3 == 3",
                                                  h4(strong("Download")),
                                                  wellPanel(
                                                    textInput("fnamexcibersort", "filename", value = "cibersort"),
                                                    downloadButton('Downloadcibersort', 'Download cibersort table'),

                                                  )),
                                 

                                 conditionalPanel("input.cPanelssession3 == 4",
                                                  wellPanel(
                                                    h4(strong("ssGSEA")),
                                                    h5("This may take some time"),
                                                    selectInput("selectssGSEAgs", "Choose gene sets",
                                                                choices=c("MsigDB.c2.cp.v6.2.symbols.gmt"),multiple = F),
                                                    actionButton("startssGSEA","start")
                                                  )),
                                 
                                 conditionalPanel("input.cPanelssession3 == 4",
                                                  h4(strong("Download")),
                                                  wellPanel(
                                                    textInput("fnamessGSEA", "filename", value = "ssGSEA"),
                                                    downloadButton('DownloadssGSEA', 'Download ssGSEA table'),
                                                    
                                                  )),
                                 
                                 
                                 conditionalPanel("input.cPanelssession3 == 5",
                                                  wellPanel(
                                                    h4(strong("Correlation")),
                                                    selectInput("selectscol1", "Choose x axis",
                                                                choices=c("A1BG"),multiple = F),
                                                    selectInput("selectscol2", "Choose y axis",
                                                                choices=c("A1CF"),multiple = F),
                                                    selectInput("selectscormethod", "Choose method",
                                                                choices=c("pearson","spearman","kendall"),
                                                                selected = "pearson",
                                                                multiple = F),
                                                  )),
                                 
                                 conditionalPanel("input.cPanelssession3 == 5",
                                                  h4(strong("Download")),
                                                  wellPanel(
                                                    textInput("fnamecorrelation", "filename", value = "correlation"),
                                                    downloadButton('Downloadcorrelation', 'Download correlation plot'),
                                                    
                                                  )),
                                 
                                 
                                 conditionalPanel("input.cPanelssession3 == 6",
                                                  wellPanel(
                                                    h4(strong("Heatmap")),
                                                    #checkboxInput("selectcluster_rows",
                                                    #              "cluster rows", TRUE),
                                                    
                                                    #checkboxInput("somevalue", "Some value", FALSE),
                                                    # selectInput("selectshow_column_names", "show column names", 
                                                    #            choices=c("YES"=T,"NO"=F),multiple = F,
                                                    #           selected = T),
                                                   # selectInput("selectshow_row_names", "show row names", 
                                                    #            choices=c("YES"=T,"NO"=F),multiple = F,
                                                    #            selected = F),
                                                    
                                                    #selectInput("selectcluster_columns", "cluster columns", 
                                                    #            choices=c("YES"=T,"NO"=F),multiple = F,
                                                     #           selected = T),
                                                    
                                                    #selectInput("selectcluster_rows", "cluster rows", 
                                                    #            choices=c("YES"=T,"NO"=F),multiple = F,
                                                    #            selected = T),
                                                    
                                                   selectInput("selectexchangecolor", "change the color of legend", 
                                                                choices=c("YES"=T,"NO"=F),multiple = F,
                                                                selected = F),
                                                   checkboxInput("selectshow_column_names", 
                                                                 "show column names", TRUE),
                                                   checkboxInput("selectshow_row_names", 
                                                                 "show row names", FALSE),
                                                   checkboxInput("selectcluster_columns", 
                                                                 "cluster columns", TRUE),
                                                   # if the argument(cluster_rows)is added ,an error will be report.
                                                   selectInput("selectupgenelist", "Upload the genelist",
                                                              choices = c("NO" = "NO","YES"="YES"),
                                                              selected="YES"
                                                    ),
                                                    conditionalPanel("input.selectupgenelist == 'YES'",
                                                                     selectInput("filegenelist",label= "choose an example or your own data", 
                                                                                 choices = c("Example"="Examplegenelist", "Your own data" = "load_my_owngenelist")),
                                                                     conditionalPanel("input.filegenelist == 'Examplegenelist'",
                                                                                      downloadButton('downloadgenelist', 'Download example')),
                                                                     conditionalPanel("input.filegenelist == 'load_my_owngenelist'",
                                                                                      fileInput('filegenelist7', 'Choose xlsx File', 
                                                                                                accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv', '.txt'))
                                                                                      
                                                                     )
                                                    ),
                                                    selectInput("selectupexprinfo", "Upload the sample info", 
                                                                choices = c("NO" = "NO","YES"="YES"),
                                                                selected="YES"
                                                    ),
                                                    conditionalPanel("input.selectupexprinfo == 'YES'",
                                                                     selectInput("fileexprinfo",label= "choose an example or your own data", 
                                                                                 choices = c("Example"="Exampleexprinfo", "Your own data" = "load_my_ownexprinfo")),
                                                                     conditionalPanel("input.fileexprinfo == 'Exampleexprinfo'",
                                                                                      downloadButton('downloadexprinfo', 'Download example')),
                                                                     conditionalPanel("input.fileexprinfo == 'load_my_ownexprinfo'",
                                                                                      fileInput('fileexprinfo8', 'Choose xlsx File', 
                                                                                                accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv', '.txt'))
                                                                                      
                                                                     )
                                                    ),
                                                
                                                    actionButton("refreshexpre", "Refresh"),
                                                    
                                                  )),
                                 
                                 conditionalPanel("input.cPanelssession3 == 6",
                                                  h4(strong("Download")),
                                                  wellPanel(
                                                    textInput("fnameheatmap", "filename", value = "heatmap"),
                                                    downloadButton('Downloadheatmap', 'Download heatmap')
                                                  )),
               
                          ),
                          
                          
                          column(10,
                                 tabsetPanel(
                                   tabPanel("Manual", htmlOutput("ReadMe4"), value =1),
                                   tabPanel("xCell", htmlOutput("pvsessionxcell"), DT::dataTableOutput("xcell",width = 1200),value = 2),
                                   tabPanel("Cibersort", htmlOutput("pvsessioncibersort"), DT::dataTableOutput("cibersort",width = 1200),value = 3),
                                   tabPanel("ssGSEA", htmlOutput("pvsessionssGSEA"), DT::dataTableOutput("ssGSEA",width = 1200),value = 4),
                                   tabPanel("Correlation", htmlOutput("pvsessionCorrelation"), plotOutput("correlation", height= 800, width = 1000), value = 5),
                                   tabPanel("Heatmap", htmlOutput("pvsessionHeatmap"), plotOutput("heatmap", height= 800, width = 1000), value = 6),
                                   
                                   id = "cPanelssession3"
                                 )
                          ),
                        )
               )
               
       
    )
)





# Define server


server <- function(input, output, session) {
  
  
  data_input1 <- reactive({
    if(input$file1 == 'Example1'){
      d2 <- read.xlsx("./example/example_mutation.xlsx")
    }
    else if(input$file1 == 'load_my_own1'){
      inFile <- input$file31
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".xlsx", inFile[1])) { d2 = read.xlsx(as.character(inFile$datapath), colNames = TRUE, rowNames = F) }
      else if(grepl(".csv", inFile[1])) { d2 = read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F, as.is = T, fill = T) }
      else if(grepl(".txt", inFile[1])) { d2 = read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, as.is = T, fill = T) }
      #openxlsx::write.xlsx(d2,"tmp.xlsx")
      }
    else 
      return(NULL)
    Dataset1 <- data.frame(d2)
    return(as.data.frame(Dataset1))
  })
  
  
  output$downloadEx1 <- downloadHandler( 
    filename <- function() {
      paste0('Example_mutation_data','.xlsx')
    },
    content <- function(file) {
      ds2 <- data_input1()
      write.xlsx(ds2, file)
    }
  )
  
  data_input2 <- reactive({
    if(input$file2 == 'Example2'){
      d2 <- read.xlsx("./example/example_info.xlsx")
    }
    else if(input$file2 == 'load_my_own2'){
      inFile <- input$file32
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".xlsx", inFile[1])) { d2 = read.xlsx(as.character(inFile$datapath), colNames = TRUE, rowNames = F) }
      else if(grepl(".csv", inFile[1])) { d2 = read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F, as.is = T, fill = T) }
      else if(grepl(".txt", inFile[1])) { d2 = read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, as.is = T, fill = T) }
    }
    else 
      return(NULL)
    Dataset2 <- data.frame(d2)
    return(as.data.frame(Dataset2))
  })
  
  
  output$downloadEx2 <- downloadHandler( 
    filename <- function() {
      paste0('Example_clinical_data','.xlsx')
    },
    content <- function(file) {
      ds2 <- data_input2()
      write.xlsx(ds2, file)
    }
  )
  
  
  data_input3 <- reactive({
    if(input$file3 == 'Example3'){
      d2 <- read.xlsx("./example/example_pathway.xlsx")
    }
    else if(input$file3 == 'load_my_own3'){
      inFile <- input$file33
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".xlsx", inFile[1])) { d2 = read.xlsx(as.character(inFile$datapath), colNames = TRUE, rowNames = F) }
      else if(grepl(".csv", inFile[1])) { d2 = read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F, as.is = T, fill = T) }
      else if(grepl(".txt", inFile[1])) { d2 = read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, as.is = T, fill = T) }
    }
    else 
      return(NULL)
    Dataset3 <- data.frame(d2)
    return(as.data.frame(Dataset3))
  })
  
  
  output$downloadEx3 <- downloadHandler( 
    filename <- function() {
      paste0('Example_pathway_data','.xlsx')
    },
    content <- function(file) {
      ds2 <- data_input3()
      write.xlsx(ds2, file)
    }
  )
  
  
  data_input4 <- reactive({
    if(input$filesur2 == 'Examplesur2'){
      d2 <- read.xlsx("./example/example_survival.xlsx")
    }
    else if(input$filesur2 == 'load_my_ownsur2'){
      inFile <- input$filesur32
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".xlsx", inFile[1])) { d2 = read.xlsx(as.character(inFile$datapath), colNames = TRUE, rowNames = F) }
      else if(grepl(".csv", inFile[1])) { d2 = read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F, as.is = T, fill = T) }
      else if(grepl(".txt", inFile[1])) { d2 = read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, as.is = T, fill = T) }
    }
    else 
      return(NULL)
    Dataset3 <- data.frame(d2)
    return(as.data.frame(Dataset3))
  })
  
  
  output$downloadExsur2 <- downloadHandler( 
    filename <- function() {
      paste0('Example_survival','.xlsx')
    },
    content <- function(file) {
      ds2 <- data_input4()
      write.xlsx(ds2, file)
    }
  )
  
  # machine learning
  data_input5 <- reactive({
    if(input$fileml2 == 'Exampleml2'){
      d2 <- read.xlsx("./example/example_mechine_learning.xlsx")
    }
    else if(input$fileml2 == 'load_my_ownml2'){
      inFile <- input$fileml32
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".xlsx", inFile[1])) { d2 = read.xlsx(as.character(inFile$datapath), colNames = TRUE, rowNames = F) }
      else if(grepl(".csv", inFile[1])) { d2 = read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F, as.is = T, fill = T) }
      else if(grepl(".txt", inFile[1])) { d2 = read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, as.is = T, fill = T) }
    }
    else 
      return(NULL)
    Dataset3 <- data.frame(d2)
    return(as.data.frame(Dataset3))
  })
  
  
  output$downloadExml2 <- downloadHandler( 
    filename <- function() {
      paste0('Example_mechine_learning','.xlsx')
    },
    content <- function(file) {
      ds2 <- data_input5()
      write.xlsx(ds2, file)
    }
  )
  
  # Expression
  data_input6 <- reactive({
    if(input$fileex2 == 'Exampleex2'){
      d2 <- read.xlsx("./example/example_expression.xlsx",rowNames = T)
    }
    else if(input$fileex2 == 'load_my_ownex2'){
      inFile <- input$fileex32
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".xlsx", inFile[1])) { d2 = read.xlsx(as.character(inFile$datapath), colNames = TRUE, rowNames = T) }
      else if(grepl(".csv", inFile[1])) { d2 = read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F, as.is = T, fill = T) }
      else if(grepl(".txt", inFile[1])) { d2 = read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, as.is = T, fill = T) }
    }
    else 
      return(NULL)
    Dataset3 <- data.frame(d2)
    return(as.data.frame(Dataset3))
  })
  
  output$downloadExex2 <- downloadHandler( 
    filename <- function() {
      paste0('Example_expression','.xlsx')
    },
    content <- function(file) {
      ds2 <- data_input6()
      write.xlsx(ds2, file,rowNames=T)
    }
  )
  
  # genelist 
  # data_input7 <- reactive({
  #   if(input$filegenelist == 'Examplegenelist'){
  #     d2 <- read.xlsx("./example/genelist.xlsx")
  #   }
  #   else if(input$filegenelist == 'load_my_owngenelist'){
  #     inFile <- input$filegenelist7
  #     if (is.null(inFile))
  #       return(NULL)
  #     else if(grepl(".xlsx", inFile[1])) { d2 = read.xlsx(as.character(inFile$datapath), colNames = TRUE, rowNames = F) }
  #     else if(grepl(".csv", inFile[1])) { d2 = read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F, as.is = T, fill = T) }
  #     else if(grepl(".txt", inFile[1])) { d2 = read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, as.is = T, fill = T) }
  #   }
  # 
  #   else 
  #     return(NULL)
  #   Dataset3 <- data.frame(d2)
  #   return(as.data.frame(Dataset3))
  # })
  data_input7 <- reactive({
    if(input$selectupgenelist=="NO"){
      return(NULL)
    }else if(input$filegenelist == 'Examplegenelist'){
      d2 <- read.xlsx("./example/genelist.xlsx")
      
    }else if(input$filegenelist == 'load_my_owngenelist'){
      inFile <- input$filegenelist7
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".xlsx", inFile[1])) { d2 = read.xlsx(as.character(inFile$datapath), colNames = TRUE, rowNames = F) }
      else if(grepl(".csv", inFile[1])) { d2 = read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F, as.is = T, fill = T) }
      else if(grepl(".txt", inFile[1])) { d2 = read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, as.is = T, fill = T) }
    }
    
    else 
      return(NULL)
    Dataset3 <- data.frame(d2)
    return(as.data.frame(Dataset3))
  })

  output$downloadgenelist <- downloadHandler( 
    filename <- function() {
      paste0('Example_expression_genelist','.xlsx')
    },
    content <- function(file) {
      ds2 <- data_input7()
      write.xlsx(ds2, file)
    }
  )
  
  # expression info
  # data_input8 <- reactive({
  #   if(input$fileexprinfo == 'Exampleexprinfo'){
  #     d2 <- read.xlsx("./example/example_expression_info.xlsx",rowNames = T)
  #   }
  #   else if(input$fileexprinfo == 'load_my_ownexprinfo'){
  #     inFile <- input$fileexprinfo8
  #     if (is.null(inFile))
  #       return(NULL)
  #     else if(grepl(".xlsx", inFile[1])) { d2 = read.xlsx(as.character(inFile$datapath), colNames = TRUE, rowNames = T) }
  #     else if(grepl(".csv", inFile[1])) { d2 = read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F, as.is = T, fill = T) }
  #     else if(grepl(".txt", inFile[1])) { d2 = read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, as.is = T, fill = T) }
  #   }
  #   else 
  #     return(NULL)
  #   Dataset3 <- data.frame(d2)
  #   return(as.data.frame(Dataset3))
  # })
  data_input8 <- reactive({
    if(input$selectupexprinfo=="NO"){
      return(NULL)
    }else if (input$fileexprinfo == 'Exampleexprinfo'){
      d2 <- read.xlsx("./example/example_expression_info.xlsx",rowNames = T)
    }
    else if(input$fileexprinfo == 'load_my_ownexprinfo'){
      inFile <- input$fileexprinfo8
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".xlsx", inFile[1])) { d2 = read.xlsx(as.character(inFile$datapath), colNames = TRUE, rowNames = T) }
      else if(grepl(".csv", inFile[1])) { d2 = read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F, as.is = T, fill = T) }
      else if(grepl(".txt", inFile[1])) { d2 = read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, as.is = T, fill = T) }
    }
    else 
      return(NULL)
    Dataset3 <- data.frame(d2)
    return(as.data.frame(Dataset3))
  })
  
  output$downloadexprinfo <- downloadHandler( 
    filename <- function() {
      paste0('Example_expression_info','.xlsx')
    },
    content <- function(file) {
      ds2 <- data_input8()
      write.xlsx(ds2, file)
    }
  )
  
  observe({
    dsnames1 <- colnames(data_input1())
    dsnames2 <- colnames(data_input2())
    dsnames4 <- colnames(data_input4())
    dsnames5 <- colnames(data_input5())
    dsnames6 <- rownames(data_input6())
    
    genelist <- unique(data_input1()[["GENE"]])
    nsample <- length(unique(data_input1()[["ORDER_ID"]]))
    vartype <-  unique(data_input1()[["VAR_TYPE"]])
    
    updateSelectInput(session, "select01", label = "Divided into two subgroups",
                      choices = c(dsnames2),
                      selected = "TMB")
    updateSelectInput(session, "select11", label = "choose clinical feature",
                      choices = c(dsnames2),
                      selected =c("GENDER","Smoker"))	
    updateSelectInput(session, "selectbox", label = "select feature",
                      choices = c(dsnames2),
                      selected = "TMB")
    updateSelectInput(session, "select_gene", label = "Gene",
                      choices =  genelist,
                      selected = genelist[1])
    updateSelectInput(session, "selectBarfeature", label = "Divided into two subgroups",
                      choices = dsnames2,
                      selected = "TMB")
    updateSelectInput(session, "select51", label = "Gene",
                      choices =  genelist,
                      selected = "TP53")
    updateSelectInput(session, "select52", label = "refseq ID",
                      choices =  dsnames1,
                      selected = "refSeqID")
    updateSelectInput(session, "selectstatfeature", label = "Select feature",
                      choices = dsnames2,
                      selected = "TMB")
    # clinical statistic
    updateSelectInput(session, "select71", label = "Divided into two subgroups",
                      choices = dsnames2,
                      selected = "TMB")
    updateSelectInput(session, "selectclifea", label = "Select clinical features",
                      choices = dsnames2,
                      selected = c("AGE","GENDER"))
    updateSliderInput(session, "numberclicutoff", label = "The value of cutoff:",
                      min=0,max =100 ,
                      value = 12, step = 1
                      )
    # statistical test 
    updateSelectInput(session, "select81", label = "Divided into two subgroups",
                      choices = dsnames2,
                      selected = "TMB")
    updateSliderInput(session, "selectmutfreq", label = "Gene mutation frequency",
                      min = 0, max =nsample, value = 5, step = 1)
    
    updateSelectInput(session, "select91", label = "Divided into two subgroups",
                      choices = dsnames2,
                      selected = "TMB")
    updateSliderInput(session, "selectmutfreq9", label = "Gene mutation frequency",
                      min = 0, max =nsample, value = 5, step = 1)
    
    # Co-Mutation
    updateSelectInput(session, "select101", label = "gene colnames",
                      choices = dsnames1,
                      selected = c("GENE"))
    updateSelectInput(session, "select102", label = "sample colnames",
                      choices = dsnames1,
                      selected = c("ORDER_ID"))
    updateSelectInput(session, "select103", label = "VAR TYPE colnames",
                      choices = dsnames1,
                      selected = c("VAR_TYPE"))
    
    # Find driver gene
    updateSelectInput(session, "select111", label = "gene colnames",
                      choices = dsnames1,
                      selected = c("GENE"))
    updateSelectInput(session, "select112", label = "sample colnames",
                      choices = dsnames1,
                      selected = c("ORDER_ID"))
    updateSelectInput(session, "select113", label = "VAR TYPE colnames",
                      choices = dsnames1,
                      selected = c("VAR_TYPE"))
    updateSelectInput(session, "select114", label = "AA CHANGE colnames",
                      choices = dsnames1,
                      selected = c("AA_CHANGE"))
    
    # Circos
    updateSelectInput(session, "select121", label = "gene colnames",
                      choices = dsnames1,
                      selected = c("GENE"))
    updateSelectInput(session, "select122", label = "var type colnames",
                      choices = dsnames1,
                      selected = c("VAR_TYPE"))
    updateSelectInput(session, "select123", label = "genomic colnames",
                      choices = dsnames1,
                      selected = c("GENOMIC"))
    updateSelectInput(session, "select124", label = "gene pair colnames",
                      choices = dsnames1,
                      selected = c("DNA_CHANGE"))
    updateSelectInput(session, "select125", label = "amp del colnames",
                      choices = dsnames1,
                      selected = c("VAR_TYPE_SX"))
    
    
    # CNV
    updateSelectInput(session, "selectVariantType", label = "Select variant type",
                      choices = vartype,
                      selected = c("CNV"))
    updateSelectInput(session, "select131", label = "gene colnames",
                      choices = dsnames1,
                      selected = c("GENE"))
    updateSelectInput(session, "select132", label = "sample colnames",
                      choices = dsnames1,
                      selected = c("ORDER_ID"))
    updateSelectInput(session, "select133", label = "cnv type colnames",
                      choices = dsnames1,
                      selected = c("DNA_CHANGE"))
    updateSelectInput(session, "select134", label = "genomic colnames",
                      choices = dsnames1,
                      selected = c("GENOMIC"))

    # survival
    updateSelectInput(session, "selecttime", label = "Time",
                      choices = dsnames4,
                      selected = "OS.time")
    updateSelectInput(session, "selectevent", label = "Event",
                      choices = dsnames4,
                      selected = "OS")
    updateSelectInput(session, "selectstrata", label = "Strata",
                      choices = dsnames4,
                      selected = "Mutation")
    # ggforest
    updateSelectInput(session, "selectforesttime", label = "Time",
                      choices = dsnames4,
                      selected = "OS.time")
    updateSelectInput(session, "selectggevent", label = "Event",
                      choices = dsnames4,
                      selected = "OS")
    
    updateSelectInput(session, "selectMultivariate", label = "Multivariate",
                      choices = dsnames4,
                      selected = c("Age", "gender","Mutation","pathologic_stage"))
    
    # mechine learning
    updateSelectInput(session, "selectLabelcol", label = "label colname",
                      choices = dsnames5,
                      selected = c("diabetes"))
    
    # correlation
    updateSelectInput(session, "selectscol1", label = "Choose x axis",
                      choices = dsnames6,
                      selected = c("A1BG"))
    updateSelectInput(session, "selectscol2", label = "Choose y axis",
                      choices = dsnames6,
                      selected = c("A1CF"))
    
  })
  
  ##profiling
  profilingforuse <- function(){
    data1 <- data_input1()
    data2 <- data_input2()
    data3 <- data_input3()
    if(input$select21=="NO"){
      res_pic <- plot_landscape(mut=data1, cli=data2, bar=input$select01,
                                bartype=input$select011,
                                feature = input$select11,
                                prefix = input$fname22,
                                n = input$number,
                                quickplot = input$selectquickplot,
                                rownamessize = input$rownamessize,
                                resam=input$selectchangecolor,
                                waterfall = input$selectwaterfall)
      
    }else{
      res_pic <- plot_pathway_cli(mut=data1, cli=data2, bar=input$select01,
                                  bartype=input$select011,
                                  feature=input$select11,
                                  prefix=input$fname22,
                                  pathway_genes = data3,
                                  nfreq=input$nfreq,
                                  rownamessize = input$rownamessize,
                                  resam = input$selectchangecolor,
                                  quickplot = input$selectquickplot,
                                  waterfall = input$selectwaterfall
                                  
      )
    }
    res_pic
    
  }
  
  
  output$profiling <- renderPlot({
    input$refresh
    isolate(profilingforuse())
  })
  
  
  output$Downloadprofiling <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fname22)
      paste(pdf_file,".pdf", sep='')
    },
    content <- function(file) {
      pdf(file, height= 10, width=12,onefile = FALSE)
      profilingforuse()
      dev.off()
      #file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    }, contentType = 'image/pdf')
  
  
  
  ##boxplot
  boxplotforuse <-  function(){
    data1 <- data_input1()
    data2 <- data_input2()
    if(input$select_type=="boxplot"){
      res_pic <- ori_boxplot(x = input$select_gene,y=input$selectbox,
                             xtype ="gene",input=data1,cli=data2)
      
    }else{
      res_pic <- ori_boxplot(x = input$select_gene,y=input$selectbox,
                             xtype ="gene",ftype="violin",input=data1,cli=data2)
    }
    print(res_pic)
    
  }
  output$boxplot <- renderPlot({
    boxplotforuse()
    
  })
  
  output$Downloadbox <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fname23)
      paste(pdf_file,'.pdf', sep='')
    },
    content <- function(file) {
      #pdf(file=paste(pdf_file,".pdf",sep="") , height= 10, width=12)
      pdf(file , height= 10, width=12,onefile = FALSE)
      boxplotforuse()
      #ggsave(filename = file)
      dev.off()
      #file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    },contentType = 'image/pdf')
  
  
  
  ##barplot
  # barplotforuse <- function(){
  #   data1 <- data_input1()
  #   if(input$selectbarplot=="NO"){
  #     if(input$selectbarplotyaxis=="percentage"){
  #       res_pic <- ori_barplot(file = data1,
  #                              id = 'ORDER_ID',
  #                              gene = 'GENE',
  #                              vartype = 'VAR_TYPE_SX',
  #                              by = NULL,
  #                              byorder = NULL,
  #                              gs = input$numberbarplot,
  #                              color = NULL,
  #                              ytype = 'percentage',
  #                              outpdf = NULL)
  #     }else if (input$selectbarplotyaxis=="counts"){
  #       res_pic <- ori_barplot(file = data1,
  #                              id = 'ORDER_ID',
  #                              gene = 'GENE',
  #                              vartype = 'VAR_TYPE_SX',
  #                              by = NULL,
  #                              byorder = NULL,
  #                              gs = input$numberbarplot,
  #                              color = NULL,
  #                              ytype = 'counts',
  #                              outpdf = NULL)
  #     }
  #     
  #   }else{
  #     mut <- data1
  #     bar <- input$selectBarfeature
  #     mut$group <- NA
  #     cutoff <- input$numberbarplotcutoff
  #     mut[[bar]] <- as.numeric(mut[[bar]])
  #     if(is.null(cutoff)){
  #       mut$group <- ifelse(mut[[bar]]>median(mut[[bar]],na.rm = T),"High","Low")  
  #     }else{
  #       mut$group <- ifelse(mut[[bar]]>cutoff,"High","Low")  
  #     }
  #     if(input$selectbarplotyaxis=="percentage"){
  #       res_pic <- ori_barplot(file = mut,
  #                              id = 'ORDER_ID',
  #                              gene = 'GENE',
  #                              vartype = 'VAR_TYPE_SX',
  #                              by = "group",
  #                              byorder = NULL,
  #                              gs = input$numberbarplot,
  #                              color = NULL,
  #                              ytype = 'percentage',
  #                              outpdf = NULL)
  #     }else if(input$selectbarplotyaxis=="counts"){
  #       res_pic <- ori_barplot(file = mut,
  #                              id = 'ORDER_ID',
  #                              gene = 'GENE',
  #                              vartype = 'VAR_TYPE_SX',
  #                              by = "group",
  #                              byorder = NULL,
  #                              gs = input$numberbarplot,
  #                              color = NULL,
  #                              ytype = 'counts',
  #                              outpdf = NULL)
  #       
  #     }
  #     
  #   }
  #   print(res_pic)
  #   
  # }

  barplotforuse <- function(){
    data1 <- data_input1()
    if(input$selectbarplot=="NO"){
      res_pic <- ori_barplot(file = data1,
                             id = 'ORDER_ID',
                             gene = 'GENE',
                             vartype = 'VAR_TYPE_SX',
                             by = NULL,
                             #byorder = NULL,
                             gs = input$numberbarplot,
                             color = NULL,
                             ytype = input$selectbarplotyaxis,
                             outpdf = NULL)
    }else{
      data2 <- data_input2()
      mut <- data1
      bar <- input$selectBarfeature
      if(input$selectBargroupty=="continuous"){
        mut$group <- NA
        cutoff <- input$numberbarplotcutoff
        mut[[bar]] <- as.numeric(mut[[bar]])
        if(is.null(cutoff)){
          mut$group <- ifelse(mut[[bar]]>median(mut[[bar]],na.rm = T),"High","Low")  
        }else{
          mut$group <- ifelse(mut[[bar]]>cutoff,"High","Low")  
        }
      }else{
        df <- unique(data2[,c("ORDER_ID",bar)])
        mut <- mut[,!colnames(mut) %in% bar]
        mut <- merge(mut,df,by="ORDER_ID",all = T)
        mut$group <- mut[[bar]]
      }
      
      res_pic <- ori_barplot(file = mut,
                             id = 'ORDER_ID',
                             gene = 'GENE',
                             vartype = 'VAR_TYPE_SX',
                             by = "group",
                             byorder = input$selectBarbyorder,
                             gs = input$numberbarplot,
                             color = NULL,
                             ytype = input$selectbarplotyaxis,
                             outpdf = NULL)
      
    }
    print(res_pic)
  }
  
 
  output$barplot <- renderPlot({
    barplotforuse()
    
  })
  
  output$Downloadbar <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fname24)
      paste(pdf_file,'.pdf', sep='')
    },
    content <- function(file) {
      pdf(file , height= 10, width=12,onefile = FALSE)
      barplotforuse()
      dev.off()
     # file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    },contentType = 'image/pdf')
  
  
  
  ##lollipop
  output$lollipop <- renderPlot({
    data1 <- data_input1()
    res_pic <- ori_lollipop(gene=input$select51,
                            mut=data1,
                            refSeqID=input$select52)
    print(res_pic)
  })
  
  output$Downloadlollipop <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fname25)
      paste(pdf_file,'.pdf', sep='')
    },
    content <- function(file) {
      pdf(file , height= 10, width=12,onefile = FALSE)
      data1 <- data_input1()
      res_pic <- ori_lollipop(gene=input$select51,
                              mut=data1,
                              refSeqID=input$select52)
      print(res_pic)
      dev.off()
      #file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    },contentType = 'image/pdf')
  
  
  
  # stat
  statforuse <- function(){
    data1 <- data_input1()
    data2 <- data_input2()
    if(input$selectstat=="NO"){
      res_pic <- CountVariants(data1,
                               id = 'ORDER_ID',
                               gene = 'GENE',
                               vartype = 'VAR_TYPE_SX',
                               varorder = c('Fusion/Rearrangement', 'Substitution/Indel', 'Gene Amplification', 'Gene Homozygous Deletion', 'Truncation')
                               # by = "group"
      )
    }else{
      bar <- input$selectstatfeature
      cli  <- data2[,c('ORDER_ID',bar),drop=F]
      mut <- merge(data1,cli,all.x = TRUE)
      mut$group <- NA
      cutoff <- input$numbercountvariantcutoff
      grouptype <- input$selectstatfeaturetype
      if(grouptype =="Continuous"){
        mut[[bar]] <- as.numeric(mut[[bar]])
        if(is.null(cutoff)){
          mut$group <- ifelse(mut[[bar]]>median(mut[[bar]],na.rm = T),"High","Low")  
        }else{
          mut$group <- ifelse(mut[[bar]]>cutoff,"High","Low")  
        }
      }else{
        mut$group <- as.character(mut[[bar]])
      }
      res_pic <- CountVariants(mut,
                               id = 'ORDER_ID',
                               gene = 'GENE',
                               vartype = 'VAR_TYPE_SX',
                               varorder = c('Fusion/Rearrangement', 'Substitution/Indel', 'Gene Amplification', 'Gene Homozygous Deletion', 'Truncation'),
                               by = "group"
      )
      return(res_pic)
    
    }
  }
  

  output$stat <- renderDataTable({
    res_pic <- statforuse()
    if(input$select61=="GeneFre"){
      res_table<- res_pic[[1]]
      res_table
    }else if(input$select61=="VarFre"){
      res_table<- res_pic[[2]]
      res_table
    }else if(input$select61=="PathwayFre"){
      res_table<- res_pic[[3]]
      res_table
    } 
  })
  
  
  output$Downloadstat <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fname26)
      paste(pdf_file,'.xlsx', sep='')
    },
    content <- function(file) {
      res_pic <- statforuse()
      write.xlsx(res_pic,file,rowNames =F)
    })
  
  #clinicalStat
  
  output$clinicalStat <- renderTable({
    data2 <- data_input2()
    re_table <- ori_statcli(input$selectclifea,input$select71,data2,
                            input$select72,
                            cutoff=input$numberclicutoff)
    re_table
  })
  
  output$Downloadclistat <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fname27)
      paste(pdf_file,'.xlsx', sep='')
    },
    content <- function(file) {
      data2 <- data_input2()
      re_table <- ori_statcli(input$selectclifea,input$select71,data2,
                              cutoff=input$numberclicutoff)
      write.xlsx(re_table,file,rowNames =F)
    }
  )
  
  # Statistical test
  
  output$statisticaltest <- renderDataTable({
    data1 <- data_input1()
    data2 <- data_input2()
    res_table <- ori_biomarker_gene(y=input$select81,
                                   mut=data1,
                                   cli= data2,
                                   n=input$selectmutfreq,
                                   ytype =input$select80,
                                   cutoff =input$numbertestcutoff)
    
    if(input$select82=="Wiltest"){
      res_table <-  res_table[[2]]
      res_table
    }else if(input$select82=="fishertest"){
      res_table <- res_table[[3]]
      res_table
      
    } 
  })
  
  output$Downloadtest <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fname28)
      paste(pdf_file,'.xlsx', sep='')
    },
    content <- function(file) {
      data1 <- data_input1()
      data2 <- data_input2()
      res_table <- ori_biomarker_gene(y=input$select81,
                                     mut=data1,
                                     cli= data2,
                                     n=input$selectmutfreq,
                                     ytype =input$select80,
                                     cutoff =input$numbertestcutoff)
      openxlsx::write.xlsx(res_table,file,rowNames =T)
    },contentType = 'text/csv')
  
  
  # fishertestplot
  output$fishertestplot <- renderPlot({
    data1 <- data_input1()
    data2 <- data_input2()
    res_table <- ori_biomarker_gene(y=input$select91,
                                    mut=data1,
                                    cli= data2,
                                    n=input$selectmutfreq9,
                                    ytype =input$select90,
                                    cutoff =input$numberfisherteplotcutoff)
    f <- res_table$fisher_test
    f$GENE <- row.names(f)
    p <- ori_fisherplot(f,OR="OR",padj ="fisher.test.Padj", GENE="GENE")
    print(p)
  })
  
  output$Downloadfisherplot <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fname29)
      paste(pdf_file,'.pdf', sep='')
    },
    content <- function(file) {
      data1 <- data_input1()
      res_table <- ori_biomarker_gene(y=input$select91,
                                      mut=data1,
                                      cli= data2,
                                      n=input$selectmutfreq9,
                                      ytype =input$select90,
                                      cutoff =input$numberfisherteplotcutoff)
      f <- res_table$fisher_test
      f$GENE <- row.names(f)
      p <- ori_fisherplot(f,OR="OR",padj ="fisher.test.Padj", GENE="GENE")
      
      pdf(file, height= 10, width=12,onefile = FALSE)
      print(p)
      dev.off()
      #file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    },contentType = 'image/pdf')
  
  # Co-Mutation
  
  output$CoMutationplot <- renderPlot({
    data1 <- data_input1()
    plot_gene_exco(mut_df=data1,gene_col=input$select101,
                   var_col=input$select103,
                   sample_col= input$select102,
                   top=input$numberCoMutation,
                   fontSize = input$select104,
                   outprefix="tmp")
  })

  output$DownloadCoMutationplot <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fname210)
      paste(pdf_file,'.pdf', sep='')
    },
    content <- function(file) {
      data1 <- data_input1()
      pdf(file, height= 6, width=8,onefile = FALSE)
      plot_gene_exco(mut_df=data1,gene_col=input$select101,
                     var_col=input$select103,
                     sample_col= input$select102,
                     top=input$numberCoMutation,
                     fontSize = input$select104,
                     outprefix="tmp")
      dev.off()
      #file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    },contentType = 'image/pdf')


  # drivergeneplot
  
  output$drivergeneplot <- renderPlot({
  data1 <- data_input1()
  drive_gene(mut_df=data1,
             gene_col=input$select111,
             sample_col=input$select112,
             AA_change_col=input$select114,
             Var_class_col=input$select113,
             outprefix="TMP",
             WriteTable = F,
             plotpdf = T)
  
})
  
  output$drivergenetable <- renderDataTable({
    data1 <- data_input1()
    res_table <- drive_gene(mut_df=data1,
                            gene_col=input$select111,
                            sample_col=input$select112,
                            AA_change_col=input$select114,
                            Var_class_col=input$select113,
                            outprefix="TMP",
                            WriteTable = T,
                            plotpdf = F)
    print(res_table)
    
  })
  
  output$DownloadDriverGeneplot <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fname211)
      paste(pdf_file,'.pdf', sep='')
    },
    content <- function(file) {
      data1 <- data_input1()
      pdf(file, height= 6, width=8,onefile = FALSE)
      drive_gene(mut_df=data1,
                 gene_col=input$select111,
                 sample_col=input$select112,
                 AA_change_col=input$select114,
                 Var_class_col=input$select113,
                 outprefix="TMP",
                 WriteTable = F,
                 plotpdf = T)
      dev.off()
      #file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    },contentType = 'image/pdf')
  
  output$DownloadDriverGenetable <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fname211)
      paste(pdf_file,'.xlsx', sep='')
    },
    content <- function(file) {
      data1 <- data_input1()
      res_table <- drive_gene(mut_df=data1,
                              gene_col=input$select111,
                              sample_col=input$select112,
                              AA_change_col=input$select114,
                              Var_class_col=input$select113,
                              outprefix="TMP",
                              WriteTable = T,
                              plotpdf = F)
      write.xlsx(res_table,file,rowNames =F)
    })
  
  #  circos
  output$circosplot <- renderPlot({
    data1 <- data_input1()
    plot_circos_gene(mut_df=data1,
                     genomic_col=input$select123,
                     mut_type_col=input$select122,
                     gene_col=input$select121,
                     outprefix="tmp",
                     plot_mut_type=input$selectMutType,
                     gene_pair_col=input$select124,
                     amp_del_col=input$select125)
  })
  
  output$DownloadCircosplot <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fname212)
      paste(pdf_file,'.pdf', sep='')
    },
    content <- function(file) {
      pdf(file , height= 8, width=8,onefile = FALSE)
      data1 <- data_input1()
      plot_circos_gene(mut_df=data1,
                       genomic_col=input$select123,
                       mut_type_col=input$select122,
                       gene_col=input$select121,
                       outprefix="tmp",
                       plot_mut_type=input$selectMutType,
                       gene_pair_col=input$select124,
                       amp_del_col=input$select125)
      dev.off()
      #file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    },contentType = 'image/pdf')
  
  #CNV
  output$CNVplot <- renderPlot({
    data1 <- data_input1()
    data2 <- data_input2()
    plot_CNV_seqment(cnv_df=data1, 
                     cli_df=data2,
                     sample_col=input$select132,
                     gene_col=input$select131,
                     genomic_col=input$select134,
                     cnv_type_col=input$select133,
                     genesize= input$select136,
                     cytobandTxtSize=input$select137,
                     outprefix="tmp",
                     top=input$select135,
                     Plotpdf = T,
                     Writetable = F
     
    )
    
  })
  
  output$CNVtable <- renderDataTable({
    data1 <- data_input1()
    data2 <- data_input2()
    plot_CNV_seqment(cnv_df=data1, 
                     cli_df=data2,
                     sample_col=input$select132,
                     gene_col=input$select131,
                     genomic_col=input$select134,
                     cnv_type_col=input$select133,
                     genesize= input$select136,
                     cytobandTxtSize=input$select137,
                     outprefix="tmp",
                     top=input$select135,
                     Plotpdf = T,
                     Writetable = T
                     
    )
    
  })
  
  output$DownloadCNVplot <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fname213)
      paste(pdf_file,'.pdf', sep='')
    },
    content <- function(file) {
      pdf(file , height= 8, width=8,onefile = FALSE)
      data1 <- data_input1()
      data2 <- data_input2()
      plot_CNV_seqment(cnv_df=data1, 
                       cli_df=data2,
                       sample_col=input$select132,
                       gene_col=input$select131,
                       genomic_col=input$select134,
                       cnv_type_col=input$select133,
                       genesize= input$select136,
                       cytobandTxtSize=input$select137,
                       outprefix="tmp",
                       top=input$select135,
                       Plotpdf = T,
                       Writetable = F
                       
      )
      dev.off()
      #file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    },contentType = 'image/pdf')
  
  output$DownloadCNVtable <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fname213)
      paste(pdf_file,'.xlsx', sep='')
    },
    content <- function(file) {
      data1 <- data_input1()
      data2 <- data_input2()
      res_table <-   plot_CNV_seqment(cnv_df=data1, 
                                      cli_df=data2,
                                      sample_col=input$select132,
                                      gene_col=input$select131,
                                      genomic_col=input$select134,
                                      cnv_type_col=input$select133,
                                      genesize= input$select136,
                                      cytobandTxtSize=input$select137,
                                      outprefix="tmp",
                                      top=input$select135,
                                      Plotpdf = T,
                                      Writetable = T
                                      
      )
      write.xlsx(res_table,file,rowNames =F)
    })

  ## TMBcompare
  output$TMBcompareplot <- renderPlot({
    data1 <- data_input1()
    tcgaCompare_ori(var = data1,
                    cohortName = input$selectcohortName,
                    medianCol = input$selectmedianCol,
                    col = c( input$selectTCGAcohort,input$selectInputcohort))
    
  })
  
  
  output$TMBcomparetable <-renderDataTable({
    data1 <- data_input1()
    res_table<- tcgaCompare_ori(var = data1,
                                cohortName = input$selectcohortName,
                                medianCol = input$selectmedianCol,
                                col = c(input$selectTCGAcohort,input$selectInputcohort))
    print(res_table)
    
  })
  
  output$DownloadTMBcompareplot <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnameTMBcompare)
      paste(pdf_file,'.pdf', sep='')
    },
    content <- function(file) {
      pdf(file , height= 10, width=8,onefile = FALSE)
      data1 <- data_input1()
      tcgaCompare_ori(var = data1,
                      cohortName = input$selectcohortName,
                      medianCol = input$selectmedianCol,
                      col = c( input$selectTCGAcohort,input$selectInputcohort))
      dev.off()
      #file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    },contentType = 'image/pdf')
  
  output$DownloadTMBcomparetable <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnameTMBcompare)
      paste(pdf_file,'.xlsx', sep='')
    },
    content <- function(file) {
      data1 <- data_input1()
      res_table <-   tcgaCompare_ori(var = data1,
                                     cohortName = input$selectcohortName,
                                     medianCol = input$selectmedianCol,
                                     col = c( input$selectTCGAcohort,input$selectInputcohort))
      write.xlsx(res_table,file,rowNames =F)
    })
  
  
  #survival
  survivalforuse <- function(){
    data4 <- data_input4()
    if(input$selectshowTable==T){
      res_pic <- ggsurvival(mx=data4, time = input$selecttime,
                            event = input$selectevent,
                            strata = input$selectstrata, 
                            showLegend = TRUE, 
                            plottitle = NULL, ytitle = 'Probability of survival', 
                            #showTable = input$selectshowTable
                            showTable = T
      )
    }else{
      res_pic <- ggsurvival(mx=data4, time = input$selecttime,
                            event = input$selectevent,
                            strata = input$selectstrata, 
                            showLegend = TRUE, 
                            plottitle = NULL, ytitle = 'Probability of survival', 
                            #showTable = input$selectshowTable
                            showTable = F
      )
    }
    print(res_pic)
  }
  
  
  output$survivalplot <- renderPlot({
    survivalforuse()
  })
  
  output$Downloadsurvival <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnamesurvival)
      paste(pdf_file,'.pdf', sep='')
    },
    content <- function(file) {
      pdf(file , height= 10, width=12,onefile = FALSE)
      survivalforuse()
      dev.off()
      #file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    },contentType = 'image/pdf')
  
  ## ggforest
  output$ggforest <- renderPlot({
    data4 <- data_input4()
    res_pic <- pggforest(mx=data4,
                         Multivar=input$selectMultivariate,
                         time = input$selectforesttime, 
                         event = input$selectggevent,
                         Stage=input$selectstage)
    print(res_pic)
    
  })
  
  output$Downloadggforest <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnamegggforest)
      paste(pdf_file,'.pdf', sep='')
    },
    content <- function(file) {
      pdf(file , height= 10, width=12,onefile = FALSE)
      data4 <- data_input4()
      res_pic <- pggforest(data4,input$selectMultivariate,
                           time = input$selectforesttime, 
                           event = input$selectggevent,
                           Stage=input$selectstage)
      print(res_pic)
      dev.off()
     # file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    },contentType = 'image/pdf')

  ## machine learning
  
  output$MLcmplot <- renderPlot({
    data5 <- data_input5()
    ori_ml(data5, 
           SplitData=input$selecteSplitData,
           Method=input$selectmethod,
           Labelcol=input$selectLabelcol,
           nCrossValidation=input$selectenCV,
           Plotimportant=F,Plotcm=T)
  })
  
  output$DownloadMLplot <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnameml)
      paste(pdf_file,'.pdf', sep='')
    },
    content <- function(file) {
      pdf(file , height= 10, width=12,onefile = FALSE)
      data5 <- data_input5()
      ori_ml(data5,SplitData=input$selecteSplitData,
             Method=input$selectmethod,
             Labelcol=input$selectLabelcol,
             nCrossValidation=input$selectenCV,
             Plotimportant=F,Plotcm=T)
      dev.off()
      # file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    },contentType = 'image/pdf')
  
  
  output$MLimportplot <- renderPlot({
    data5 <- data_input5()
    ori_ml(data5,SplitData=input$selecteSplitData,
           Method=input$selectmethod,
           Labelcol=input$selectLabelcol,
           nCrossValidation=input$selectenCV,
           Plotimportant=T,Plotcm=F)
  })
  
  output$DownloadMLimportplot <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnameml)
      paste(pdf_file,'_importance.pdf', sep='')
    },
    content <- function(file) {
      pdf(file , height= 10, width=12,onefile = FALSE)
      data5 <- data_input5()
      ori_ml(data5,SplitData=input$selecteSplitData,
             Method=input$selectmethod,
             Labelcol=input$selectLabelcol,
             nCrossValidation=input$selectenCV,
             Plotimportant=T,Plotcm=F)
      dev.off()
      # file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    },contentType = 'image/pdf')
  
  #xCell
  
  # observeEvent(input$startxCell, {
  #   session$sendCustomMessage(type = 'testmessage',
  #                             message = 'Please wait,it takes some time')
  # })
  
  observeEvent(input$startxCell, {
    showNotification(paste("Please wait,it takes some time"), duration = 5)
  })
  
  res_table_xcell <- eventReactive(input$startxCell, {
    data6 <- data_input6()
    xCellAnalysis(data6)
  })

  output$xcell <- renderDataTable({
      res_table <- res_table_xcell()
      ##print(res_table)
    })
  
  output$Downloadxcell <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnamexcell)
      paste(pdf_file,'.xlsx', sep='')
    },
    content <- function(file) {
      res_table <- res_table_xcell()
      write.xlsx(res_table,file,rowNames =T)
    }
  )
 
  # xcellresult <- reactive({
  #   data6 <- data_input6()
  #   res_table<-xCellAnalysis(data6)
  #   return(res_table)
  # 
  # })
  # 
  # output$xcell <- renderDataTable({
  #   res_table <- xcellresult()
  #   #print(res_table)
  # })
  # 
  # output$Downloadxcell <- downloadHandler(
  #   filename <- function() {
  #     pdf_file <<- as.character(input$fnamexcell)
  #     paste(pdf_file,'.xlsx', sep='')
  #   },
  #   content <- function(file) {
  #     res_table <- xcellresult()
  #     write.xlsx(res_table,file,rowNames =T)
  #   }
  # )

  
  # output$xcell <- renderDataTable({
  #   data6 <- data_input6()
  #   res_table<-xCellAnalysis(data6)
  #   #print(res_table)
  # })
  # 
  # output$Downloadxcell <- downloadHandler(
  #   filename <- function() {
  #     pdf_file <<- as.character(input$fnamexcell)
  #     paste(pdf_file,'.xlsx', sep='')
  #   },
  #   content <- function(file) {
  #     data6 <- data_input6()
  #     res_table<- xCellAnalysis(data6)
  #     write.xlsx(res_table,file,rowNames =T)
  #   }
  # )
  
  #cibersort
  observeEvent(input$startcibersort, {
    showNotification(paste("Please wait,it takes some time"), duration = 5)
  })
  
  res_table_ciber <- eventReactive(input$startcibersort, {
    data6 <- data_input6()
    data6$name <- rownames(data6)
    data6 <- data6[,c(ncol(data6),1:(ncol(data6)-1))]
    res_table <- CIBERSORT("./data/LM22.txt",data6,perm = 100, absolute = F)
  })
  
  output$cibersort <- renderDataTable({
    res_table <- res_table_ciber()
    #print(res_table)
  })
  
  output$Downloadcibersort <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnamexcibersort)
      paste(pdf_file,'.xlsx', sep='')
    },
    content <- function(file) {
      res_table <- res_table_ciber()
      write.xlsx(res_table,file,rowNames =T)
    }
  )
  
  # ssGSEA
  observeEvent(input$startssGSEA, {
    showNotification(paste("Please wait,it takes some time"), duration = 5)
  })
  
  res_table_ssGSEA <- eventReactive(input$startssGSEA, {
    data6 <- data_input6()
    data61 <- as.matrix(data6)
    gene_sets = fgsea::gmtPathways("./data/MsigDB.c2.cp.v6.2.symbols.gmt")
    res_table<-GSVA::gsva(expr = data61, gset.idx.list = gene_sets, min.sz=10,max.sz=500,method = 'ssgsea', 
                          ssgsea.norm=F,
                          verbose = TRUE)
  })
  
  output$ssGSEA <- renderDataTable({
    res_table<-res_table_ssGSEA()
    #print(res_table)
  })
  
  output$DownloadssGSEA <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnamessGSEA)
      paste(pdf_file,'.xlsx', sep='')
    },
    content <- function(file) {
      res_table<-res_table_ssGSEA()
      write.xlsx(res_table,file,rowNames =T)
    }
  )
  
  # correlation
  output$correlation <- renderPlot({
    data6 <- data_input6()
    data6.1 <- as.data.frame(t(data6))
    res_pic <- ori_scattercor_plot(data6.1,
                                   x=input$selectscol1,
                                   y=input$selectscol2,
                                   method=input$selectscormethod)
    print(res_pic)
  })
  
  output$Downloadcorrelation <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnamecorrelation)
      paste(pdf_file,'.pdf', sep='')
    },
    content <- function(file) {
      pdf(file , height= 10, width=12,onefile = FALSE)
      data6 <- data_input6()
      data6.1 <- as.data.frame(t(data6))
      res_pic <- ori_scattercor_plot(data6.1,
                                     x=input$selectscol1,
                                     y=input$selectscol2,
                                     method=input$selectscormethod)
      print(res_pic)
      dev.off()
      # file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    },contentType = 'image/pdf')
  
  # Heatmap
  Heatmapforuse <- function(){
    expr <- data_input6() 
    genelist <- data_input7()
    info <- data_input8()
    p <- ori_heatmap(expr=expr,
                     cli = info,
                     genelist = genelist,
                     resam=input$selectexchangecolor,
                     show_column_names=input$selectshow_column_names,
                     show_row_names=input$selectshow_row_names,
                     cluster_columns=input$selectcluster_columns,
                     #cluster_rows=input$selectcluster_rows,
                     
    )
    print(p)
    
  }
  
  output$heatmap <- renderPlot({
    input$refreshexpre
    isolate(Heatmapforuse())
    
  })
  
  output$Downloadheatmap <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnameheatmap)
      paste(pdf_file,'.pdf', sep='')
    },
    content <- function(file) {
      pdf(file , height= 10, width=12,onefile = FALSE)
      Heatmapforuse()
      dev.off()
      # file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    },contentType = 'image/pdf')
  

  ## ReadMe
  output$ReadMe1 <- renderUI({
    str00 <- paste("&emsp;")
    str0 <- paste("示例")
    str1 <- paste("&emsp; 1.在公司内部网络环境下或者在公司外部使用公司配置的笔记本电脑，打开浏览器，输入 http://cotest.origimed.com/Himalaya/，即可使用分析工具")
    str2 <- paste("&emsp; 2.示例文件包括 突变文件，临床信息，如果需要展示profiling在信号通路上的突变，则需要信号通路文件")
    str3 <- paste("&emsp; 3.在示例文件下，点击不同模块，调整对应参数,即可展示结果，点击下载按钮，输入文件名，
                  即可将结果保存在目的文件夹下，若不输入文件名，则保存为默认文件")
    
    str21 <- paste("数据格式")
    str22 <- paste("&emsp; 1.建议上传文件格式为.xlsx 或者下载示例数据，将文件内容替换，再重新上传，注意文件解密")
    str23 <- paste("&emsp; 2.上传的任何文件都不允许使用中文字符。样本ID列名必须是 ORDER_ID，突变基因列记录的列名是GENE，突变类型列记录列名是VAR_TYPE_SX")
    str24 <- paste("&emsp; 3.上传的文件，若存在与示例文件列内容相同，但列名不同时，建议将列名更改为与示例文件列名相同，尤其注意列名大小写问题")

    str31 <- paste("分析模块")
    str32 <- paste("&emsp; 1. Profiling &emsp;需要突变文件和临床文件，若需要展示在信号通路上的突变，则需要信号通路文件，
                   该部分可以根据临床信息展示在profiling上方，若要改变默认颜色，选择YES，并且点击 Refresh可以改变颜色，
                   该部分选择参数发生改变时，都需要点击refresh。  
                   若没有临床文件，允许仅使用突变文件，注意：在上传突变文件后，
                   在Input your clinical file 选项，选择Your own data。若上传文件样本数过多，可选择 quick plot ，然后调整参数，出图。
                   建议一般情况下quick plot 选择为NO
                   ")
    str33 <- paste("&emsp; 2. Boxplot &emsp;可以根据基因突变状态进行分组，并展示两组之间TMB或者其他列信息，
                   注意这里要求展示的变量是连续性变量，展示方式可以是boxplot或者小提琴图")
    str34 <- paste("&emsp; 3. Barplot &emsp; 可以展示前n个高频突变基因，并且可以选择cutoff进行分组展示。
                   y轴可以进行调整为percentage或者counts")
    str35 <- paste("&emsp; 4. Lollipop &emsp;选择对应基因的Lollipop图")
    str36 <- paste("&emsp; 5. Count variant &emsp;统计基因的突变情况，包括是基因突变频率，基因不同突变类型突变频率，在信号通路上突变频率，
                   下载时，将三个表格一同下载。")
    str37 <- paste("&emsp; 6. Clinical statistics &emsp;根据分组统计对应临床信息")
    str38 <- paste("&emsp; 7. Statistics test &emsp;包括wilcoxon test和 t test 以及 fisher test和chisq test，
                   可以分别查看，下载时是一同下载")
    str39 <- paste("&emsp; 8. Fishertest plot &emsp;可以选择使用分组的信息，并且选择cutoff，以及筛选高频突变基因，
                   Gene mutation frequency是指突变样本个数")
    str310 <- paste("&emsp; 9. Co-mutation plot &emsp;可以选择高频突变基因，默认是前30个高频突变基因。可以调整gene name的大小。
                   Gene mutation frequency是指突变样本个数")
    str311 <- paste("&emsp;10. Find driver gene &emsp; 基于突变位点的位置信息，使用oncodriveCLUST算法识别driver gene，可分别下载plot和table")
    str312 <- paste("&emsp;11. Circos &emsp; 可以选择在Circos上展示的mutation type，可以单独选择一种或两种，在这里默认的是展示snv，cnv，fusion三种")
    str313 <- paste("&emsp;12. CNV &emsp; 展示CNV Amp和Del在染色体上突变情况和突变高频基因，图和表可分别下载。")
    str314 <- paste("&emsp;13. TMBcompare &emsp;与TCGA cohort的TMB进行比较，可以调整cohort在图中展示的名字以及颜色，图和表可分别下载")
    
    
    HTML(paste(str00,h5(strong(str0)), str1, str2, str3,str00,h5(strong(str21)),str22,str23,str24,str00,h5(strong(str31)),
         str32,str33,str34,str35,str36,str37,str38,str39,str310,str311,str312,str313,str314,sep = '<br/>'))
  })
  
  output$ReadMe2 <- renderUI({
    str00 <- paste("&emsp;")
    str0 <- paste("示例")
    str1 <- paste("&emsp; 1.在公司内部网络环境下，打开浏览器，输入 http://10.10.174.10:3838/Himalaya/，则可使用分析工具")
    str2 <- paste("&emsp; 2.示例文件包括 生存分析文件，即包括生存时间，生存状态，分组信息")
    str3 <- paste("&emsp; 3.在示例文件下，点击不同模块，调整对应参数,即可展示结果，点击下载按钮，输入文件名，
                  即可将结果保存在目的文件夹下，若不输入文件名，则保存为默认文件")
    
    str21 <- paste("数据格式")
    str22 <- paste("&emsp; 1.建议上传文件格式为.xlsx 或者下载示例数据，将文件内容替换，再重新上传，注意文件解密")
    str23 <- paste("&emsp; 2.上传的任何文件都不允许使用中文字符。")
    str24 <- paste("&emsp; 3.上传的文件，若存在与示例文件列内容相同，但列名不同时，建议将列名更改为与示例文件列名相同，尤其注意列名大小写问题")
    
    str31 <- paste("分析模块")
    str32 <- paste("&emsp; 1.Survival &emsp;需要三列信息，即生存时间，生存状态，分组信息。可选择是否出现统计表格")
    str33 <- paste("&emsp; 2.ggforest &emsp;需要三列信息，即生存时间，生存状态，以及选择进行分析的因素，若forest plot中HR过大，则可调整分析的因素")
    
    HTML(paste(str00,h5(strong(str0)), str1, str2, str3,str00,h5(strong(str21)),str22,str23,str24,str00,h5(strong(str31)),
               str32,str33,sep = '<br/>'))
  })

  output$ReadMe3 <- renderUI({
    str00 <- paste("&emsp;")
    str0 <- paste("示例")
    str1 <- paste("&emsp; 1.在公司内部网络环境下，打开浏览器，输入 http://10.10.174.10:3838/Himalaya/，则可使用分析工具")
    str2 <- paste("&emsp; 2.示例文件包括特征列和标签列")
    str3 <- paste("&emsp; 3.在示例文件下，点击Machine Learning，调整对应参数,即可展示结果，点击下载按钮，输入文件名，
                  即可将结果保存在目的文件夹下，若不输入文件名，则保存为默认文件")
    
    str21 <- paste("数据格式")
    str22 <- paste("&emsp; 1.建议上传文件格式为.xlsx 或者下载示例数据，将文件内容替换，再重新上传，注意文件解密")
    str23 <- paste("&emsp; 2.上传的任何文件都不允许使用中文字符。")
    str24 <- paste("&emsp; 3.上传的文件，若存在与示例文件列内容相同，但列名不同时，建议将列名更改为与示例文件列名相同，尤其注意列名大小写问题")
    
    str31 <- paste("分析模块")
    str32 <- paste("&emsp; MechineLearning &emsp;需要标签列信息，即分组信息。在训练模型时，可选择使用的不同的算法，
                    包括Random Forest（随机森林），Support Vector Machine（支持向量机），Logistic Regression(逻辑回归)，Naive Bayes（朴素贝叶斯），
                    k-Nearest  Neighbors（K近邻），Decision Tree（决策树），默认是Random Forest。可选择交叉验证倍数，默认为5倍交叉验证。
                    可选择将数据切割成训练集和验证集的比例，这里默认是将全部数据的0.8划分为训练集。
                   ")
    HTML(paste(str00,h5(strong(str0)), str1, str2, str3,str00,h5(strong(str21)),str22,str23,str24,str00,h5(strong(str31)),
               str32,sep = '<br/>'))
  })
  
  
  output$ReadMe4 <- renderUI({
    str00 <- paste("&emsp;")
    str0 <- paste("示例")
    str1 <- paste("&emsp; 1.在公司内部网络环境下，打开浏览器，输入 http://10.10.174.10:3838/Himalaya/，则可使用分析工具")
    str2 <- paste("&emsp; 2.示例文件为表达文件 行是基因，列是样本。")
    str3 <- paste("&emsp; 3.在示例文件下，点击不同模块，调整对应参数,即可展示结果，点击下载按钮，输入文件名，
                  即可将结果保存在目的文件夹下，若不输入文件名，则保存为默认文件")
    str4 <- paste("&emsp; 4.注意：xCell，Cibersort，ssGSEA需要运行一定时间，在运行过程中，耐心等待运行结果，请勿来回调整界面或者重复点击start")
    
    str21 <- paste("数据格式")
    str22 <- paste("&emsp; 1.建议上传文件格式为.xlsx 或者下载示例数据，将文件内容替换，再重新上传，注意文件解密")
    str23 <- paste("&emsp; 2.上传的任何文件都不允许使用中文字符。")
    str24 <- paste("&emsp; 3.上传的文件，若存在与示例文件列内容相同，但列名不同时，建议将列名更改为与示例文件列名相同，尤其注意列名大小写问题")
    
    str31 <- paste("分析模块")
    str32 <- paste("&emsp; 1.xCell &emsp;导入表达文件，点击start按钮，即可开始进行xcell分析，在这里是默认的gene signature是64种免疫细胞，分析需要运行一定的时间")
    str33 <- paste("&emsp; 2.Cibersort &emsp;导入表达文件，点击start按钮，即可开始进行cibersort分析，在这里是默认的gene signature是22种免疫细胞，分析需要运行一定的时间")
    str34 <- paste("&emsp; 3.ssGSEA &emsp;导入表达文件，点击start按钮，即可开始进行ssGSEA分析，在这里是默认的gene sets是MsiDB.c2.cp.v6.2，分析需要运行一定的时间")
    str35 <- paste("&emsp; 4.Correlation &emsp;导入表达文件，可选择任意两列计算相关计算，计算方法包括pearson，spearman，kendall。")
    str36 <- paste("&emsp; 5.Heatmap &emsp;可上传需要展示的genelist和sample info，如无genelist文件，默认展示前50个基因。若要改变默认颜色，选择YES，并且点击 Refresh改变颜色。
                   该部分选择参数发生改变时，都需要点击refresh。")
    HTML(paste(str00,h5(strong(str0)), str1, str2, str3, str4,str00,h5(strong(str21)),str22,str23,str24,str00,h5(strong(str31)),
               str32,str33,str34,str35,str36,sep = '<br/>'))
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
