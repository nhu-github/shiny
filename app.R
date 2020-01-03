#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# add change test 
library(shiny)
library(openxlsx)
library(ggplot2)
library(shinythemes)
library(DT)

#source("./script/ori_profiling_shiny.R")
source("./script/ori_profiling_shiny_quickplot.R")
source("./script/ori_boxplot_shiny.R")
source("./script/ori_barplot.R")
source("./script/CountVariants.R")
source("./script/ori_tools.R")
source("./script/lollipop.R")
source("./script/ori_lollipop_shiny.R")
source("./script/ori_statcli_shiny.R")
source("./script/ori_fisherplot_shiny.R")
source("./script/ori_biomarker_gene_shiny.R")
source("./script/ori_survival.R")
source("./script/gene_ex_co_shiny.r")
source("./script/drive_gene_shiny.r")
source("./script/circus_shiny.r")
source("./script/CNV_shiny.r")

# Define UI for application 
ui <- fluidPage(
    theme = shinytheme('cerulean'),
    navbarPage("Pipeline", 
               
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
                                                        selectInput("select11", "Choose clinical feature", 
                                                                    choices=c("GENDER"),multiple = T),
                                                        
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
                                                        
    
                                                        sliderInput("number", 
                                                                    label = "The number of first top mutant genes:",
                                                                    min = 10, max = 100, value = 30, step = 1
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
                                                                                     choices = c("TMB"))
                                                        ),
                                                        sliderInput("numberbarplotcutoff", 
                                                                    label = "The value of cutoff:",
                                                                    min = 5, max = 50, value = 12, step = 1
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
  
  
  
  observe({
    dsnames1 <- colnames(data_input1())
    dsnames2 <- colnames(data_input2())
    dsnames4 <- colnames(data_input4())
    
    genelist <- unique(data_input1()[["GENE"]])
    nsample <- length(unique(data_input1()[["ORDER_ID"]]))
    vartype <-  unique(data_input1()[["VAR_TYPE"]])
    
    updateSelectInput(session, "select01", label = "Divided into two subgroups",
                      choices = c(dsnames1),
                      selected = "TMB")
    updateSelectInput(session, "select11", label = "choose clinical feature",
                      choices = c(dsnames2),
                      selected =c("GENDER","Smoker"))	
    updateSelectInput(session, "selectbox", label = "select feature",
                      choices = c(dsnames1),
                      selected = "TMB")
    updateSelectInput(session, "select_gene", label = "Gene",
                      choices =  genelist,
                      selected = genelist[1])
    updateSelectInput(session, "selectBarfeature", label = "Divided into two subgroups",
                      choices = dsnames1,
                      selected = "TMB")
    updateSelectInput(session, "select51", label = "Gene",
                      choices =  genelist,
                      selected = "TP53")
    updateSelectInput(session, "select52", label = "refseq ID",
                      choices =  dsnames1,
                      selected = "refSeqID")
    updateSelectInput(session, "selectstatfeature", label = "Divided into two subgroups",
                      choices = dsnames1,
                      selected = "TMB")
    updateSelectInput(session, "select71", label = "Divided into two subgroups",
                      choices = dsnames2,
                      selected = "TMB")
    updateSelectInput(session, "selectclifea", label = "Select clinical features",
                      choices = dsnames2,
                      selected = c("AGE","GENDER"))
    
    # statistical test 
    updateSelectInput(session, "select81", label = "Divided into two subgroups",
                      choices = dsnames1,
                      selected = "TMB")
    updateSliderInput(session, "selectmutfreq", label = "Gene mutation frequency",
                      min = 0, max =nsample, value = 5, step = 1)
    
    updateSelectInput(session, "select91", label = "Divided into two subgroups",
                      choices = dsnames1,
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
    
    
  })
  
  ##profiling
  profilingforuse <- function(){
    data1 <- data_input1()
    data2 <- data_input2()
    data3 <- data_input3()
    if(input$select21=="NO"){
      res_pic <- plot_landscape(mut=data1, cli=data2, bar=input$select01, 
                                feature = input$select11,
                                prefix = input$fname22,
                                n = input$number,
                                quickplot = input$selectquickplot,
                                resam=input$selectchangecolor,
                                waterfall = input$selectwaterfall)
      
    }else{
      res_pic <- plot_pathway_cli(mut=data1, cli=data2, bar=input$select01, 
                                  feature=input$select11,
                                  prefix=input$fname22,
                                  pathway_genes = data3,
                                  nfreq=0.03,
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
    if(input$select_type=="boxplot"){
      res_pic <- ori_boxplot(x = input$select_gene,y=input$selectbox,
                             xtype ="gene",input=data1)
      
    }else{
      res_pic <- ori_boxplot(x = input$select_gene,y=input$selectbox,
                             xtype ="gene",ftype="violin",input=data1)
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
  barplotforuse <- function(){
    data1 <- data_input1()
    if(input$selectbarplot=="NO"){
      if(input$selectbarplotyaxis=="percentage"){
        res_pic <- ori_barplot(file = data1,
                               id = 'ORDER_ID',
                               gene = 'GENE',
                               vartype = 'VAR_TYPE_SX',
                               by = NULL,
                               byorder = NULL,
                               gs = input$numberbarplot,
                               color = NULL,
                               ytype = 'percentage',
                               outpdf = NULL)
      }else if (input$selectbarplotyaxis=="counts"){
        res_pic <- ori_barplot(file = data1,
                               id = 'ORDER_ID',
                               gene = 'GENE',
                               vartype = 'VAR_TYPE_SX',
                               by = NULL,
                               byorder = NULL,
                               gs = input$numberbarplot,
                               color = NULL,
                               ytype = 'counts',
                               outpdf = NULL)
      }
      
    }else{
      mut <- data1
      bar <- input$selectBarfeature
      mut$group <- NA
      cutoff <- input$numberbarplotcutoff
      mut[[bar]] <- as.numeric(mut[[bar]])
      if(is.null(cutoff)){
        mut$group <- ifelse(mut[[bar]]>median(mut[[bar]],na.rm = T),"High","Low")  
      }else{
        mut$group <- ifelse(mut[[bar]]>cutoff,"High","Low")  
      }
      if(input$selectbarplotyaxis=="percentage"){
        res_pic <- ori_barplot(file = mut,
                               id = 'ORDER_ID',
                               gene = 'GENE',
                               vartype = 'VAR_TYPE_SX',
                               by = "group",
                               byorder = NULL,
                               gs = input$numberbarplot,
                               color = NULL,
                               ytype = 'percentage',
                               outpdf = NULL)
      }else if(input$selectbarplotyaxis=="counts"){
        res_pic <- ori_barplot(file = mut,
                               id = 'ORDER_ID',
                               gene = 'GENE',
                               vartype = 'VAR_TYPE_SX',
                               by = "group",
                               byorder = NULL,
                               gs = input$numberbarplot,
                               color = NULL,
                               ytype = 'counts',
                               outpdf = NULL)
        
      }
      
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
  output$stat <- renderDataTable({
    data1 <- data_input1()
    if(input$selectstat=="NO"){
      res_pic <- CountVariants(data1,
                               id = 'ORDER_ID',
                               gene = 'GENE',
                               vartype = 'VAR_TYPE_SX',
                               varorder = c('Fusion/Rearrangement', 'Substitution/Indel', 'Gene Amplification', 'Gene Homozygous Deletion', 'Truncation')
                               # by = "group"
      )
    }else{
      mut <- data1
      bar <- input$selectstatfeature
      mut$group <- NA
      cutoff <- NULL
      mut[[bar]] <- as.numeric(mut[[bar]])
      if(is.null(cutoff)){
        mut$group <- ifelse(mut[[bar]]>median(mut[[bar]],na.rm = T),"High","Low")  
      }else{
        mut$group <- ifelse(mut[[bar]]>cutoff,"High","Low")  
      }
      res_pic <- CountVariants(mut,
                               id = 'ORDER_ID',
                               gene = 'GENE',
                               vartype = 'VAR_TYPE_SX',
                               varorder = c('Fusion/Rearrangement', 'Substitution/Indel', 'Gene Amplification', 'Gene Homozygous Deletion', 'Truncation'),
                               by = "group"
      )
    }
    
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
      data1 <- data_input1()
      if(input$selectstat=="NO"){
        res_pic <- CountVariants(data1,
                                 id = 'ORDER_ID',
                                 gene = 'GENE',
                                 vartype = 'VAR_TYPE_SX',
                                 varorder = c('Fusion/Rearrangement', 'Substitution/Indel', 'Gene Amplification', 'Gene Homozygous Deletion', 'Truncation')
                                 # by = "group"
        )
      }else{
        mut <- data1
        bar <- input$selectstatfeature
        mut$group <- NA
        cutoff <- input$numbercountvariantcutoff
        mut[[bar]] <- as.numeric(mut[[bar]])
        if(is.null(cutoff)){
          mut$group <- ifelse(mut[[bar]]>median(mut[[bar]],na.rm = T),"High","Low")  
        }else{
          mut$group <- ifelse(mut[[bar]]>cutoff,"High","Low")  
        }
        res_pic <- CountVariants(mut,
                                 id = 'ORDER_ID',
                                 gene = 'GENE',
                                 vartype = 'VAR_TYPE_SX',
                                 varorder = c('Fusion/Rearrangement', 'Substitution/Indel', 'Gene Amplification', 'Gene Homozygous Deletion', 'Truncation'),
                                 by = "group"
        )
      }
      write.xlsx(res_pic,file,rowNames =F)
    })
  
  #clinicalStat
  
  output$clinicalStat <- renderTable({
    data2 <- data_input2()
    re_table <- ori_statcli(input$selectclifea,input$select71,data2,
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
    res_table <- ori_biomarker_gene(input$select81,data1,input$selectmutfreq,
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
      res_table <- ori_biomarker_gene(input$select81,data1,input$selectmutfreq,
                                      cutoff =input$numbertestcutoff)
      openxlsx::write.xlsx(res_table,file,rowNames =T)
    },contentType = 'text/csv')
  
  
  # fishertestplot
  output$fishertestplot <- renderPlot({
    data1 <- data_input1()
    res_table <- ori_biomarker_gene(input$select91,data1,input$selectmutfreq9,
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
      res_table <- ori_biomarker_gene(input$select91,data1,input$selectmutfreq9,
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
    data2 <- data1[data1$VAR_TYPE==input$selectVariantType,]
    plot_CNV_seqment(cnv_df=data2, 
                     sample_col=input$select132,
                     gene_col=input$select131,
                     genomic_col=input$select134,
                     cnv_type_col=input$select133,
                     outprefix="tmp",
                     top=10,
                     Plotpdf = T,
                     Writetable = F
    )
    
  })
  
  output$CNVtable <- renderDataTable({
    data1 <- data_input1()
    data2 <- data1[data1$VAR_TYPE==input$selectVariantType,]
    
    # res_table <- plot_CNV_seqment(cnv_df=data2, 
    #                               sample_col=input$select132,
    #                               gene_col=input$select131,
    #                               genomic_col=input$select134,
    #                               cnv_type_col=input$select133,
    #                               outprefix="tmp",
    #                               top=10,
    #                               Plotpdf = T,
    #                               Writetable = T
    # )
    # print(res_table)
    plot_CNV_seqment(cnv_df=data2, 
                     sample_col=input$select132,
                     gene_col=input$select131,
                     genomic_col=input$select134,
                     cnv_type_col=input$select133,
                     outprefix="tmp",
                     top=10,
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
      data2 <- data1[data1$VAR_TYPE==input$selectVariantType,]
      plot_CNV_seqment(cnv_df=data2, 
                       sample_col=input$select132,
                       gene_col=input$select131,
                       genomic_col=input$select134,
                       cnv_type_col=input$select133,
                       outprefix="tmp",
                       top=10,
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
      data2 <- data1[data1$VAR_TYPE==input$selectVariantType,]
      res_table <-  plot_CNV_seqment(cnv_df=data2, 
                                     sample_col=input$select132,
                                     gene_col=input$select131,
                                     genomic_col=input$select134,
                                     cnv_type_col=input$select133,
                                     outprefix="tmp",
                                     top=10,
                                     Plotpdf = T,
                                     Writetable = T
      )
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
  
  
  ## ReadMe
  output$ReadMe1 <- renderUI({
    str00 <- paste("&emsp;")
    str0 <- paste("示例")
    str1 <- paste("&emsp; 1.在公司内部网络环境下，打开浏览器，输入 http://10.10.174.10:3838/pipeline1/，即可使用分析工�?")
    str2 <- paste("&emsp; 2.示例文件包括 突变文件，临床信息，如果需要展示profiling在信号通路上的突变，则需要信号通路文件")
    str3 <- paste("&emsp; 3.在示例文件下，点击不同模块，调整对应参数,即可展示结果，点击下载按钮，输入文件名，
                  即可将结果保存在目的文件夹下，若不输入文件名，则保存为默认文�?")
    
    str21 <- paste("数据格式")
    str22 <- paste("&emsp; 1.建议上传文件格式�?.xlsx 或者下载示例数据，将文件内容替换，再重新上传，注意文件解密")
    str23 <- paste("&emsp; 2.上传的任何文件都不允许使用中文字符。样本ID列名必须�? ORDER_ID，突变基因列记录的列名是GENE，突变类型列记录列名是VAR_TYPE_SX")
    str24 <- paste("&emsp; 3.上传的文件，若存在与示例文件列内容相同，但列名不同时，建议将列名更改为与示例文件列名相同，尤其注意列名大小写问题")

    str31 <- paste("分析模块")
    str32 <- paste("&emsp; 1. Profiling &emsp;需要突变文件和临床文件，若需要展示在信号通路上的突变，则需要信号通路文件�?
                   该部分可以根据临床信息展示在profiling上方，若要改变默认颜色，选择YES，并且点�? Refresh可以改变颜色�?
                   该部分选择参数发生改变时，都需要点击refresh。若没有临床文件，允许仅使用突变文件，注意：在上传突变文件后�?
                   在Input your clinical file 选项，选择Your own data。若上传文件样本数过多，可选择 quick plot ，然后调整参数，出图�?
                   建议一般情况下quick plot 选择为NO
                   ")
    str33 <- paste("&emsp; 2. Boxplot &emsp;可以根据基因突变状态进行分组，并展示两组之间TMB或者其他列信息�?
                   注意这里要求展示的变量是连续性变量，展示方式可以是boxplot或者小提琴�?")
    str34 <- paste("&emsp; 3. Barplot &emsp; 可以展示前n个高频突变基因，并且可以选择cutoff进行分组展示�?
                   y轴可以进行调整为percentage或者counts")
    str35 <- paste("&emsp; 4. Lollipop &emsp;选择对应基因的Lollipop�?")
    str36 <- paste("&emsp; 5. Count variant &emsp;统计基因的突变情况，包括是基因突变频率，基因不同突变类型突变频率，在信号通路上突变频率，
                   下载时，将三个表格一同下载�?")
    str37 <- paste("&emsp; 6. Clinical statistics &emsp;根据分组统计对应临床信息")
    str38 <- paste("&emsp; 7. Statistics test &emsp;包括wilcoxon test�? t test以及 fisher test和chisq test�?
                   可以分别查看，下载时是一同下�?")
    str39 <- paste("&emsp; 8. Fishertest plot &emsp;可以选择使用分组的信息，并且选择cutoff，以及筛选高频突变基因，
                   Gene mutation frequency是指突变样本个数")

    HTML(paste(str00,h5(strong(str0)), str1, str2, str3,str00,h5(strong(str21)),str22,str23,str24,str00,h5(strong(str31)),
         str32,str33,str34,str35,str36,str37,str38,str39, sep = '<br/>'))
  })
  
  output$ReadMe2 <- renderUI({
    str00 <- paste("&emsp;")
    str0 <- paste("示例")
    str1 <- paste("&emsp; 1.在公司内部网络环境下，打开浏览器，输入 http://10.10.174.10:3838/pipeline1/，则可使用分析工�?")
    str2 <- paste("&emsp; 2.示例文件包括 生存分析文件，即包括生存时间，生存状态，分组信息")
    str3 <- paste("&emsp; 3.在示例文件下，点击不同模块，调整对应参数,即可展示结果，点击下载按钮，输入文件名，
                  即可将结果保存在目的文件夹下，若不输入文件名，则保存为默认文�?")
    
    str21 <- paste("数据格式")
    str22 <- paste("&emsp; 1.建议上传文件格式�?.xlsx 或者下载示例数据，将文件内容替换，再重新上传，注意文件解密")
    str23 <- paste("&emsp; 2.上传的任何文件都不允许使用中文字符�?")
    str24 <- paste("&emsp; 3.上传的文件，若存在与示例文件列内容相同，但列名不同时，建议将列名更改为与示例文件列名相同，尤其注意列名大小写问题")
    
    str31 <- paste("分析模块")
    str32 <- paste("&emsp; 1.Survival &emsp;需要三列信息，即生存时间，生存状态，分组信息。可选择是否出现统计表格")
    str33 <- paste("&emsp; 2.ggforest &emsp;需要三列信息，即生存时间，生存状态，以及选择进行分析的因素，若forest plot中HR过大，则可调整分析的因素")

    HTML(paste(str00,h5(strong(str0)), str1, str2, str3,str00,h5(strong(str21)),str22,str23,str24,str00,h5(strong(str31)),
               str32,str33,sep = '<br/>'))
  })
  
  
  
  

}


# Run the application 
shinyApp(ui = ui, server = server)
