library(shiny)
library(networkD3)
library(shinybusy)
library(tictoc)


# Define UI for data upload app ----
ui <- fluidPage(  
  tabsetPanel(
    tabPanel("Construction of dataset",
             # App title ----
             titlePanel("Uploading File"),
             
             # Sidebar layout with input and output definitions ----
             sidebarLayout(
               
               # Sidebar panel for inputs ----
               sidebarPanel(
                 
                 # Input: Select a file ----
                 fileInput("file1", "Choose CSV File",
                           multiple = TRUE,
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                 
                 # Horizontal line ----
                 tags$hr(),
                 
                 # Input: Checkbox if file has header ----
                 checkboxInput("header", "Header", TRUE),
                 
                 # Input: Select separator ----
                 radioButtons("sep", "Separator",
                              choices = c(Comma = ",",
                                          Semicolon = ";",
                                          Tab = "\t"),
                              selected = ","),
                 
                 # Input: Select quotes ----
                 radioButtons("quote", "Quote",
                              choices = c(None = "",
                                          "Double Quote" = '"',
                                          "Single Quote" = "'"),
                              selected = '"'),
                 
                 # Horizontal line ----
                 tags$hr(),
                 
                 # Input: Select number of rows to display ----
                 radioButtons("disp", "Display",
                              choices = c(Head = "head",
                                          All = "all"),
                              selected = "head")
                 
               ),
               
               # Main panel for displaying outputs ----
               mainPanel(checkboxGroupInput("variables","Choose the variables you want to enrich your dataset",
                                            choiceNames = list("Terms(Go.Db)","Orthology(Kegg.db)","Motif(Kegg.db)","Pathway(Kegg.db)","Band(Biomart)","Transcript count(Biomart)","Gene biotype(Biomart)"),
                                            choiceValues = list("terms","orthology","motif","pathway","band","transcriptcount","genebiotype")
                                            
               ),
               
               textOutput("checkbox"),
               add_busy_spinner(spin = "fading-circle"),
               
               # Output: Data file ----
               textOutput("data1"),
               DT::dataTableOutput("clusters"),
               tags$hr(),
               
               textOutput("data2"),
               
               DT::dataTableOutput("refinal"),
               tags$hr(),
               
               textOutput("data3"),
               
               DT::dataTableOutput("re2")
               )
               
               
             )
    ),
    tabPanel("Parameters",
             titlePanel("Select the parametres of your analysis"),
             sidebarLayout(
               sidebarPanel( 
                 radioButtons(inputId = "algorithmcenter", label = "Choose if you want an optimal number of clusters or you want to decide the number of clusters", choices = c("optimal","mynumber"),selected =  "optimal") ,
                 sliderInput("numalgorithm", "Select the number of custers", min=2, max=15,value = 2, step=1) ,
                 radioButtons(inputId = "algorithm",label = "Select the algorithm you want to use(note that this choice continues to all the further analysis)",choices = c("kmeans","kmodes","kproto"), selected = "kmeans")
                 #textOutput("checkmaxnumber")
               ),
               mainPanel(
                 add_busy_spinner(spin = "fading-circle"),
                 
                 tags$hr(),
                 
                 textOutput("kmeanss"), 
                 #flowLayout(
                 
                 
                 tableOutput("clust1"),
                 #plotOutput("kmeanssil"),
                 #),
                 
                 tags$hr(),
                 textOutput("kmodess"),
                 #flowLayout(
                 tableOutput("clust2"),
                 #plotOutput("kmodesdif")
                 #),
                 
                 tags$hr(),
                 textOutput("kprotoss"),
                 
                 #flowLayout(
                 tableOutput("clust3"),
                 plotOutput("kprotodif")
                 #)
                 
                 
                 
               )
             )
    ),
    tabPanel("Categorical data analysis",selectInput("barplotvar", label = "Select the variable you want to have the Barplot of a scecific cluster",choices = c("terms","orthology","motif","pathway1","pathway2","pathway3","pathway4","pathway5","pathwaynum","band","transcript_count","gene_biotype"), selected = NULL) ,
             selectInput("numcharact", label = "Choose the number of the most frequent characteristics for the variable you chose",choices = c(2:10)),
             selectInput("numclust", label = "Choose the number of the cluster you want to see the analysis ",choices = c(1:30),selected = 1),
             
             plotOutput("barplotclust1"),
             selectInput("columnbox", label = "Choose the column of arthmetic variable you want to analyze by reference to the categorical variable",choices = c(2:140),selected = 2),
             
             plotOutput("boxplot1"),
             textOutput("namesofvar"),
             DT::dataTableOutput("refrequent"),
             
             # plotOutput("barplotclust2"),
             # plotOutput("barplotclust3"),
             # plotOutput("barplotclust4"),
             # plotOutput("barplotclust5"),
             # plotOutput("barplotclust6"),
             # plotOutput("barplotclust7"),
             # plotOutput("barplotclust8"),
             # plotOutput("barplotclust9"),
             
             # plotOutput("barplotclust10"),
             # plotOutput("barplotclust11"),
             #plotOutput("barplotclust12"),
             #plotOutput("barplotclust13"),
             #plotOutput("barplotclust14"),
             # plotOutput("barplotclust15"),
             tableOutput("recateg")
             
    ),
    
    tabPanel("Arithmetical data analysis",selectInput("columnchoice", label = "Select a column of numericaldata 1(select only numerical variables)",choices = c(2:130)) ,
             selectInput("columnchoice1", label = "Select a column of numericaldata 2(select only numerical variables)",choices = c(3:130)) ,
             
             plotOutput("clustofnumer"),
             plotOutput("boxplot2"),
             plotOutput("colourofclust")
             # tableOutput("check")
    ),
    
    tabPanel("Upset diagram",titlePanel("The most frequent values of two variables that are common for the genes of the cluster(you chosed in categorical data analysis)  of the clustering") ,
             selectInput("barplotvar123", label = "Select the first categorical variable for the upset diagram",choices = c("terms","orthology","motif","pathway1","pathway2","pathway3","pathway4","pathway5","pathwaynum","band","transcript_count","gene_biotype"),selected = "terms"),
             selectInput("barplotvar1234", label = "Select the second categorical variable for the upset diagram",choices = c("terms","orthology","motif","pathway1","pathway2","pathway3","pathway4","pathway5","pathwaynum","band","transcript_count","gene_biotype"),selected = "orthology"),
             plotOutput("upset"),
             textOutput("colnames1"),
             textOutput("colnames2"),
             textOutput("colnames3"),
             textOutput("colnames4"),
             textOutput("colnames5"),
             textOutput("colnames6")
             
    ),
    
    
    tabPanel("Sankey Diagram",titlePanel("The amount of the genes that belong in the same cluster for every clustering"),
             #textInput("ensemblid", "Put an ansembl id", "value"),
             sankeyNetworkOutput("sankeydiagram"),
             DT::dataTableOutput("resankeysou"),
             DT::dataTableOutput("resankeynames")),
    
    tabPanel("Entropy", 
             selectInput("barplotvar12", label = "Select the variable you want to see the barpolt for the entropy of every cluster",choices = c("terms","orthology","motif","pathway1","pathway2","pathway3","pathway4","pathway5","pathwaynum","band","transcript_count","gene_biotype"), selected = NULL) ,
             plotOutput("entropybarplot"),
             textOutput("entropytext"),
             textOutput("refrequentnames")
             
    ),
    tabPanel("Cohesion", 
             #tableOutput("clust4"),
             plotOutput("plotgower"),
             textOutput("colorss")
             
    )
    
  )
)
