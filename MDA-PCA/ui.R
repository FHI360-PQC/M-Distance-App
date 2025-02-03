library(shinydashboard)
library(DT)
library(mdatools)
library(data.table)

#to modify app
#library(rsconnect)
#deployApp()

#sidebar design
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Reference Data", tabName = "reference", icon = icon("table")),
    menuItem("Sample Data", tabName = "sample", icon = icon("chart-line"))
  )
)

#_____________________________________________________________________________________


#body design
body <- dashboardBody(
  tabItems(

    tabItem(tabName = "reference",
            tabsetPanel(
              tabPanel("Upload",
                       fluidRow(
                         column(3,
                                box(width = 12, title = "Upload Reference Data",
                                    # Input: File type ----
                                    selectInput("filetype1", "File type",
                                                 choices = c(LabSpec = "labspec",
                                                             Tellspec = "tellspec",
                                                             #CalFile = "calfile",
                                                             TI_App = "tiapp",
                                                             Compressed = 'compressed'
                                                 ),
                                                 selected = "tellspec"),

                                    # Input: Select a file ----
                                    fileInput("file1", "Choose File",
                                              multiple = TRUE,
                                              accept = c("text/csv",
                                                         "text/comma-separated-values,text/plain",
                                                         ".csv")),

                                    # Add background file
                                      checkboxInput("bg", tags$b("Background Correct")),
                                      conditionalPanel(
                                        condition = "input.bg == '1'",
                                        fileInput("background", "Background File",
                                                  multiple = TRUE,
                                                  accept = c("text/csv",
                                                             "text/comma-separated-values,text/plain",
                                                             ".csv"))
                                      ),
                                        
                                      checkboxInput("logtrans", tags$b("Log Transform")),
                                    
                                    selectInput("trunc_ref", "Truncation",
                                                choices = c("None","Manual"),
                                                selected = "None"
                                    ),
                                    conditionalPanel(
                                      condition = "input.trunc_ref == 'Manual'",
                                      numericInput("start_nm_ref","Min Wavelength",950),
                                      numericInput("end_nm_ref", "Max Wavelength",1650)
                                    ),
                                    
                                    # Horizontal line ----
                                    tags$hr()

                                )#end box

                         ), #end column 1
                         column(8,
                                tags$style(type="text/css",
                                           ".shiny-output-error { visibility: hidden; }",
                                           ".shiny-output-error:before { visibility: hidden; }"
                                ),

                                box(width=14,
                                    plotOutput("spectra"),
                                    downloadButton("downloadPlot", "Plot"),
                                    downloadButton("downloadData", "Data")
                                ),

                         ) #end column 2
                       ), #end fluidRow

                       fluidRow(
                         box(
                           column(12, align='center',
                                  DTOutput("contents"), style = "overflow-x: scroll;"
                           ),
                           width="100%"
                         )
                       ) # end fluidRow
              ), #end tabPanel

              tabPanel("PCA",
                       fluidRow(
                         shinyjs::useShinyjs(),
                         column(3,
                                box(width = 12, title = "PCA Parameters",
                                    checkboxInput("center", tags$b("Center Data"), value = TRUE),
                                    checkboxInput("scale", tags$b("Scale Data"), value = FALSE),
                                    checkboxInput("residuals", tags$b("Use Residuals"), value = TRUE),
                                    numericInput("components", tags$b("Principal Components"), value = 5),
                                    #downloadButton("downloadCalFile", "Cal File")
                                )#end box

                         ), #end column 1

                         column(3,
                                verbatimTextOutput("pcaSummary")
                         ),

                         column(5,
                                plotOutput("screeplot"),
                                br(),
                                plotOutput("residuals")
                         )

                       ) #end fluidRow
              ) #end tabPanel
              ,

              tabPanel("M-Distance",
                       fluidRow(
                         box(
                           title="Quartile Summary of M-Dists by Outcome",
                           DTOutput("mdist_sum1"), style = "overflow-x: scroll;"
                         ),
                         column(5,
                                box(width = "50%", title = "M-Dist Options",
                                    numericInput("mdistCI", tags$b("Confidence Interval"), value = 0.95, max = 0.99, step = .01),
                                    downloadButton("data1Download", "Download Analysis")
                                )
                         ) #end column 1
                       ), #end fluidRow

                       fluidRow(
                         box(
                           title="Means per Sample",
                           column(12, align='center',
                                  DTOutput("mdist"), style = "overflow-x: scroll;"
                           ),
                           width="100%"
                         )
                       ), # end fluidRow

                       fluidRow(
                         box(
                           title = "Means Per Column",
                           selectInput('mdist_iden1', tags$b("Column to Average by:"),
                                       choices=c('Outcome'), selected='Outcome'),
                           column(12, align='center',
                                  DTOutput("mdist_grp1"), style = "overflow-x: scroll;"
                           ),
                           width = "100%"
                         )
                       ), # end fluidRow
              ), #end tabPanel
              
              tabPanel("Normality",
                       fluidRow(
                         box (title='M-Distance Normality',
                              selectInput('ref.normalityplottype',
                                          'Plot Type:',
                                          choices=c('PDF','CDF'),
                                          selected='PDF'),
                              selectInput('ref.normalityselect',
                                          'Plot Normality by:',
                                          choices=c('Upload Reference Data'),
                                          selected='Upload Reference Data'),
                              plotOutput("ref.normalityplot"),
                              downloadButton("downloadRefNormPlot", "Plot"),
                              width='100%'
                         ) # end box
                       ) # end fluidRow
              ) # end tabPanel
            ) #end tabsetPanel
    ), #end tabItem

    #_____________________________________________________________________________________

    #body design

    tabItem(tabName = "sample",
            tabsetPanel(
              tabPanel("Upload",
                       fluidRow(
                         column(3,
                                box(width = 12, title = "Upload Sample Data",

                                    # Input: File type ----
                                    selectInput("filetype2", "File type",
                                                 choices = c(LabSpec = "labspec",
                                                             Tellspec = "tellspec",
                                                             TI_App = "tiapp",
                                                             Compressed = 'compressed'
                                                 ),
                                                 selected = "tellspec"),

                                    # Input: Select a file ----
                                    fileInput("file2", "Choose CSV File",
                                              multiple = TRUE,
                                              accept = c("text/csv",
                                                         "text/comma-separated-values,text/plain",
                                                         ".csv")),
                                    checkboxInput("bg2", tags$b("Background Correct")),
                                    
                                    # Tellspec Background / Log Transform
                                    conditionalPanel(
                                      condition = "input.bg2 == '1'",
                                      fileInput("background2", "Background File",
                                                multiple = TRUE,
                                                accept = c("text/csv",
                                                           "text/comma-separated-values,text/plain",
                                                           ".csv"))
                                    ),
                                    checkboxInput("logtrans2", tags$b("Log Transform")),
                                    
                                    selectInput("trunc_samp", "Truncation",
                                                choices = c("None","Auto","Manual"),
                                                selected = "Auto"
                                    ),
                                    conditionalPanel(
                                      condition = "input.trunc_samp == 'Manual'",
                                      numericInput("start_nm_samp","Min Wavelength",950),
                                      numericInput("end_nm_samp", "Max Wavelength",1650)
                                    ),
                                   
                                )#end box

                         ), #end column 1

                         column(8,
                                tags$style(type="text/css",
                                           ".shiny-output-error { visibility: hidden; }",
                                           ".shiny-output-error:before { visibility: hidden; }"
                                ),

                                box(width=14,
                                    plotOutput("spectra2"),
                                    downloadButton("downloadPlot2", "Plot"),
                                    downloadButton("downloadData2", "Data")
                                ),

                         ) #end column 2
                       ), #end fluidRow

                       fluidRow(
                         box(
                           column(12, align='center',
                                  DTOutput("contents2"), style = "overflow-x: scroll;"
                           ),
                           width="100%"
                         )
                       ) # end fluidRow

              ), #end tabPanel

              tabPanel("M-Distance",
                       fluidRow(
                         box(
                           title="Quartile Summary of M-Dists by Outcome",
                           DTOutput("mdist_sum2"), style = "overflow-x: scroll;"
                         ),
                         column(5,
                                box(width = "50%", title = "M-Dist Options",
                                    downloadButton("data2Download", "Download Analysis")
                                )#end box
                         ) #end column 1
                       ), #end fluidRow

                       fluidRow(
                         box(
                           title="Means per Sample",
                           column(12, align='center',
                                  DTOutput("mdist2"), style = "overflow-x: scroll;"
                           ),
                           width="100%"
                         )
                       ), # end fluidRow

                       fluidRow(
                         box(
                           title = "Means Per Column",
                           selectInput('mdist_iden2', tags$b("Column to Average by:"),
                                       choices=c('Outcome'), selected='Outcome'),
                           column(12, align='center',
                                  DTOutput("mdist_grp2"), style = "overflow-x: scroll;"
                           ),
                           width="100%"
                         )
                       ), # end fluidRow

              ), #end tabPanel
              
              tabPanel("Normality",
                       fluidRow(
                          box (title='M-Distance Normality',
                               selectInput('sample.normalityplottype',
                                           'Plot Type:',
                                           choices=c('PDF','CDF'),
                                           selected='PDF'),
                               selectInput('sample.normalityselect',
                                           'Plot Normality by:',
                                           choices=c('Upload Sample Data'),
                                           selected='Upload Sample Data'),
                               plotOutput("sample.normalityplot"),
                               downloadButton("downloadSampNormPlot", "Plot"),
                               width='100%'
                               ) # end box
                       ) # end fluidRow
              ) # end tabPanel
            ) #end tabsetPanel
    ) #end tabItem
  ) #end tabItems
) #end dasbhoardBody


#_____________________________________________________________________________________

#App arguments
dashboardPage(
  dashboardHeader(title = "PCA MDR"),
  sidebar,
  body
)
