library(shiny)
library(dplyr)
library(readr)
library(DT)
library(chemometrics)
library(mdatools)
library(factoextra)
library(ggplot2)
library(ggdist)

server <- function(input, output, session) {

  #increase file size uploads
  options(shiny.maxRequestSize=30*2048^2)

  #stop app when closing
  session$onSessionEnded(function() {
    stopApp()
  })

  observe({
    if(input$filetype1=="calfile"){
      shinyjs::disable("center")
      shinyjs::disable("scale")
      shinyjs::disable("residuals")
      shinyjs::disable("components")
      shinyjs::disable("mdistCI")
    } else {
      shinyjs::enable("center")
      shinyjs::enable("scale")
      shinyjs::enable("residuals")
      shinyjs::enable("components")
      shinyjs::enable("mdistCI")
    }
  })
  
  # Metadata colname selectInput Observes
  observeEvent(input$file1, { #Reference selectInput Options
    col.options <- colnames(getMdist()[[1]])
    updateSelectInput(session, 'mdist_iden1', "Column to Average by:", 
                      choices = col.options, selected='Outcome')
    updateSelectInput(session, "ref.normalityselect", "Plot Normality by:", 
                      choices = col.options, selected='Outcome')
  })
  observeEvent(input$file2, { #Sample selectInput Options
    col.options <- colnames(getMdist2()[[1]])
    updateSelectInput(session, 'mdist_iden2', "Column to Average by:", 
                      choices = col.options, selected='Outcome')
    updateSelectInput(session, "sample.normalityselect", "Plot Normality by:", 
                      choices = col.options, selected='Outcome')
  })

  #_________________________________________________________________________________________________________________________________
  #REFERENCE DATA

  ## M-Distance Calculation
  getMdist <- reactive({
    A <- getRefData()[[1]]
    
    ctr <- if(input$center == TRUE) {
      "TRUE"
    } else {
      "FALSE"
    }
    
    scl <- if(input$scale == TRUE) {
      "TRUE"
    } else {
      "FALSE"
    }
    
    #scale data
    A <- prep.autoscale(A, center = as.logical(ctr), scale =  as.logical(scl))
    
    # #savitzky-golay filtering
    # if(input$filter == TRUE) {
    #   A <- prep.savgol(A, width = as.numeric(input$window), porder = as.numeric(input$porder), dorder = as.numeric(input$sgolay))
    # }
    
    #PCA
    pca <- pca(A, ncomp = input$components)
    PCs <- pca(A)
    
    #get scores
    S <- pca$calres$scores
    
    #get factors
    F <- pca$loadings
    
    #sum of squared residuals
    R <- pca$calres$Q[,input$components]
    
    #mean center error
    Rc <- scale(R, center = as.logical(ctr), scale = as.logical(scl))
    
    #append error
    Sr <- if(input$residuals == TRUE) {
      cbind(S[,1:input$components], Rc)
    } else {
      S[,1:input$components]
    }
    
    #mahalanobis distance
    M <- (t(Sr)%*%Sr)/(dim(A)[1]-1)
    
    #inverse mahalanobis matrix
    M_1 <- solve(M, tol = 1e-70)
    
    #squared m-distance
    D2 <- Sr %*% M_1 %*% t(Sr)
    
    #m-distance
    D <- diag(D2)^(1/2)
    
    #mahalanobis w/RMSG
    RMSG <- ((sum(D^2)/(dim(A)[1]-1)))^(1/2)
    mdist2 <- D/RMSG
    
    output <- as.data.frame(cbind("M-Dist" = mdist2, "Residual" = R))
    
    #M-distance cutoff calculation
    n <- dim(A)[1]
    k <- input$components
    
    f_crit <- qf(input$mdistCI, df1=k, df2=(n-k))
    Dm <- k*(n-1)
    Dma <- f_crit*Dm
    Dmax <- Dma/(n-k)
    DM <- sqrt(Dmax)
    
    cutoffDiff <- mdist2-DM
    PF <- ifelse(mdist2<DM, "Pass", "Fail")
    mean_mdist <- mean(mdist2)
    
    output <- as.data.frame(cbind("Sample"=row.names(output), "M-Dist" = mdist2, "Residual" = R, "Factors" = k, "Cutoff" = DM,
                                  "Difference" = cutoffDiff, "Outcome" = PF))
    
    identifiers <- getRefData()[[2]]
    output <- left_join(output,identifiers, by="Sample")
    
    #Objects returned
    return(list(output,A, pca, df, Rc, M_1, RMSG, PCs, DM, k, mean_mdist))
    
  }) #end getMdist reactive
  
  ## Upload Tab
  getRefData <- reactive({
    dataset <- load_all_filetypes(input$filetype1, input$file1$datapath, input$file1$name)
    data_cols <- grepl("^(?=.)([+-]?([0-9]*)(\\.([0-9]+))?)$", colnames(dataset), perl=TRUE)
    
    data <- dataset[data_cols]
    metadata <- dataset[!data_cols]
    
    # Average Together duplicate SampleNames
    data <- cbind(metadata[1], data)
    colnames(data)[1] <- "SampleName"
    data <- aggregate(.~SampleName, data, FUN=mean, na.action = NULL, na.rm = TRUE)
    
    # Make data have sample names as rownames
    rownames(data) <- data[,1]
    data <- data[,2:ncol(data)]
    
    # Truncate data if desired
    if (input$trunc_ref == "Manual") {
      data <- truncate_by_wavelength(data, input$start_nm_ref, input$end_nm_ref)
    }
    
    # Collapse metadata into same number of rows as data
    colnames(metadata)[1] <- "Sample"

    if (ncol(metadata) > 1) {metadata <- aggregate(.~Sample, metadata, FUN=first, na.action = NULL)}
    else{metadata <- unique(metadata)}
    
    # Apply Background and/or log transform
    if(input$bg == TRUE) {
      
      # Get Background Data
      bg <- load_all_filetypes(input$filetype1, input$background$datapath, input$background$name)
      bg <- bg[data_cols]
      
      # Truncate bg if desired
      if (input$trunc_ref == "Manual") {
        bg <- truncate_by_wavelength(bg, input$start_nm_ref, input$end_nm_ref)
      }
      
      bgmeans <- colMeans(bg)
      
      for(i in 1:length(bgmeans)) data[,i] <- data[,i]/bgmeans[i]
    }
    
    if(input$logtrans == TRUE) {
      data <- log10(1/data)
    }
    
    return(list(data,metadata))
  })
  
  #data table output
  output$contents <- renderDT({
    df <- getRefData()[[1]]
    return(df)
    options = list(scrollX = TRUE)
  })
  
  #spectra output
  output$spectra <- renderPlot(spectraPlot())
  
  spectraPlot <- function()({
    df <- getRefData()[[1]]
    A <- df
    mdaplot(A, type = 'l')
  })
  
  #data download handler (spectral data)
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("ref_spectral_data.csv")
    },
    content = function(file) {
      df <- getRefData()[[1]]
      scan_sample <- rownames(df)
      df <- cbind(scan_sample, df)
      write.csv(df, file, row.names = FALSE)
    }
  ) #end data download handler
  
  #plot download handler (spectral plot)
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0("ref_spectral_plot.png")
    },
    content = function(file) {
      png(file, width=1000)
      spectraPlot()
      dev.off()
    }
  ) #end plot download handler
  
  
  ## PCA Tab
  
  #pca output
  output$pcaSummary <- renderPrint({
    summary(getMdist()[[8]])
  }) #end output
  
  #screeplot output
  output$screeplot <- renderPlot({
    plotVariance(getMdist()[[3]], type = 'b', show.labels = TRUE, labels = 'values')
    
  })
  
  #screeplot output
  output$residuals <- renderPlot({
    plotResiduals(getMdist()[[3]], show.labels = TRUE)
  })
  
  #data download handler (calfile data)
  #output$downloadCalFile <- downloadHandler(
  #  filename = function() {
  #    paste0("ref_calfile.csv")
  #  },
  #  content = function(file) {
  #    df <- getRefData()[[1]]
  #    scan_sample <- rownames(df)
  #    df <- cbind(scan_sample, df)
  #    
  #    output <- as.data.frame(cbind(df, "Factors" = input$components, "Center" = input$center,
  #                                  "Scale" = input$scale, "Residuals" = input$residuals, "CI" = input$mdistCI))
  #    write.csv(output, file, row.names = FALSE)
  #  }
  #) #end data download handler
  
  
  ## M-Dist Tab
  
  #mdist output
  output$mdist <- renderDataTable(datatable(getMdist()[[1]]) %>%
                                    formatStyle("Outcome", target="row",
                                                backgroundColor = styleEqual(
                                                  levels = c("Pass","Fail"),
                                                  values = c("#e6ffee","#ffe6e6"))
                                    ) %>% formatStyle("Outcome", fontWeight = "bold")
                                  )
  
  output$mdist_sum1 <- renderDT({summarize_mdist_table(getMdist()[[1]])})
  
  output$mdist_grp1 <- renderDataTable(datatable(group_mdist_table(getMdist()[[1]],
                                                                   input$mdist_iden1)) %>%
                                         formatStyle("Outcome", fontWeight = "bold") %>%
                                         formatStyle("Outcome", target="row",
                                                     backgroundColor = styleEqual(
                                                       levels = c("Pass","Fail"),
                                                       values = c("#e6ffee","#ffe6e6")))
  )
  
  #data download handler (m-dist)
  output$data1Download <- downloadHandler(
    filename = "Reference.csv",
    content = function(file) {
      mean_mdist <- getMdist()[[11]]
      date_reported <- as.character(Sys.Date())
      head <- rbind("Date Reported"=date_reported,"Mean M-Distance"=mean_mdist,"")
      write.table(head, file,sep=",",quote=FALSE,row.names = TRUE,col.names = FALSE)
      write.table(getMdist()[[1]],file,append=TRUE,sep=",", row.names = FALSE)
      write.table("",file,append=TRUE,sep=",",row.names = FALSE, col.names = FALSE)
      write.table("Quantiles per Outcome",file,append=TRUE,sep=",",
                  row.names = FALSE, col.names = FALSE)
      write.table(as.matrix(t(c("","Passes", "Fails"))),
                  file, append=TRUE,sep=",",row.names = FALSE, col.names = FALSE)
      write.table(summarize_mdist_table(getMdist()[[1]]),file,append=TRUE,sep=",", row.names = TRUE, col.names = FALSE)
      write.table("",file,append=TRUE,sep=",",row.names = FALSE, col.names = FALSE)
      write.table(paste0("Per ", input$mdist_iden1, " Data"),file,append=TRUE,sep=",", row.names = FALSE, col.names = FALSE)
      write.table(group_mdist_table(getMdist()[[1]], input$mdist_iden1),file,append=TRUE,sep=",", row.names = FALSE)
    }
  ) #end data download handler
  
  # Reference Normality Plot
  output$ref.normalityplot <- renderPlot(NormalityPlot(getMdist()[[1]], input$ref.normalityselect,
                                                       kind=input$ref.normalityplottype))
  
  #plot download handler (normality plot)
  output$downloadRefNormPlot <- downloadHandler(
    filename = function() {
      paste0("ref_normality_plot.png")
    },
    content = function(file) {
      ggsave(file, NormalityPlot(getMdist()[[1]], input$ref.normalityselect,
                                 kind=input$ref.normalityplottype), 
             width=2500, height=1125, units='px')
    }
  ) #end plot download handler
  
  
  #_________________________________________________________________________________________________________________________________
  #SAMPLE DATA
  
  ## M-Dist Calculation
  
  #Analyze sample data
  getMdist2 <- reactive({
    
    #get reference data frame
    df <- getMdist()[[2]]
    
    #get existing PCA model
    pca <- getMdist()[[3]]
    
    # Get center/scale
    ctr <- if(input$center == TRUE) {
      "TRUE"
    } else {
      "FALSE"
    }
    
    scl <- if(input$scale == TRUE) {
      "TRUE"
    } else {
      "FALSE"
    }
    
    #get scores
    S <- pca$calres$scores
    
    #get factors
    F <- pca$loadings
    
    #sum of squared residuals
    R <- pca$calres$Q[,input$components]
    Rc <- scale(R, center = TRUE, scale = FALSE)
    
    #get new data
    newdata <- getSampleData()[[1]]
    
    A <- newdata
    
    #scale new data
    center <- attr(df,"prep:center")
    scale <- attr(df,"prep:scale")
    
    A2 <- scale(A, center, scale) %*% pca$loadings
    
    
    #get error
    E <- scale(A,center, scale) - A2%*%t(pca$loadings)
    
    
    #subtract mean of training group residuals
    E2 <- rowSums(E^2)
    E3 <- E2-attr(Rc,"scale")
    
    #append error to new spectral data
    A3 <- cbind(A2,E3)
    
    #get mahalanobis matrix
    M_1 <- getMdist()[[6]]
    
    #
    D2 <- A3 %*% M_1 %*% t(A3)
    
    #
    D <- diag(D2)^(1/2)
    
    #get RMSG
    RMSG <- getMdist()[[7]]
    
    #mahalanobis w/RMSG
    mdist2 <- D/RMSG
    
    DM <- getMdist()[[9]]
    
    k <- getMdist()[[10]]
    
    cutoffDiff <- mdist2-DM
    PF <- ifelse(mdist2<DM, "Pass", "Fail")
    
    mean_mdist <- mean(mdist2)
    
    output <- as.data.frame(cbind("Sample"=row.names(A3), "M-Dist" = mdist2, "Residual" = E2, "Factors" = k, "Cutoff" = DM,
                                  "Difference" = cutoffDiff, "Outcome" = PF))
    
    identifiers <- getSampleData()[[2]]
    output <- left_join(output,identifiers, by="Sample")
    
    return(list(output,A, pca, mean_mdist))
    
  }) #end getMdist2 reactive
  
  ## Upload Tab
  getSampleData <- reactive({
    dataset <- load_all_filetypes(input$filetype2, input$file2$datapath, input$file2$name)
    data_cols <- grepl("^(?=.)([+-]?([0-9]*)(\\.([0-9]+))?)$", colnames(dataset), perl=TRUE)
    
    data <- dataset[data_cols]
    metadata <- dataset[!data_cols]
    
    # Average Together duplicate SampleNames
    data <- cbind(metadata[1], data)
    colnames(data)[1] <- "SampleName"
    data <- aggregate(.~SampleName, data, FUN=mean, na.action = NULL, na.rm = TRUE)
    
    # Make data have sample names as rownames
    rownames(data) <- data[,1]
    data <- data[,2:ncol(data)]
    
    # Truncate data if desired
    if (input$trunc_samp == "Manual") {
      data <- truncate_by_wavelength(data, input$start_nm_samp, input$end_nm_samp)
    }
    else if (input$trunc_samp == "Auto" & !is.null(input$file1)) {
      data <- auto_truncate(getRefData()[[1]], data)
    }
    
    # Collapse metadata into same number of rows as data
    colnames(metadata)[1] <- "Sample"
    
    if (ncol(metadata) > 1) {metadata <- aggregate(.~Sample, metadata, FUN=first, na.action = NULL)}
    else{metadata <- unique(metadata)}
    
    # Apply Background and/or log transform
    if(input$bg2 == TRUE) {
      
      # Get Background Data
      bg <- load_all_filetypes(input$filetype2, input$background2$datapath, input$background2$name)
      bg <- bg[data_cols]
      
      # Truncate bg if desired
      if (input$trunc_samp == "Manual") {
        bg <- truncate_by_wavelength(bg, input$start_nm_samp, input$end_nm_samp)
      }
      else if (input$trunc_samp == "Auto" & !is.null(input$file1)) {
        bg <- auto_truncate(getRefData()[[1]], bg)
      }
      
      bgmeans <- colMeans(bg)
      
      for(i in 1:length(bgmeans)) data[,i] <- data[,i]/bgmeans[i]
    }
    
    if(input$logtrans2 == TRUE) {
      data <- log10(1/data)
    }
    
    return(list(data,metadata))
  })
  
  #data table output
  output$contents2 <- renderDT({
    df <- getSampleData()[[1]]
    return(df)
    options = list(scrollX = TRUE)
  })
  
  #spectra output
  output$spectra2 <- renderPlot(spectraPlot2())
  
  spectraPlot2 <- function()({
    A <- getSampleData()[[1]]
    #A <- df
    mdaplot(A, type = 'l')
  })
  
  
  
  ## M-Dist Tab
  
  #mdist output
  output$mdist2 <- renderDataTable(datatable(getMdist2()[[1]]) %>%
                                     formatStyle("Outcome", target="row",
                                                 backgroundColor = styleEqual(
                                                   levels = c("Pass","Fail"),
                                                   values = c("#e6ffee","#ffe6e6"))
                                     ) %>% formatStyle("Outcome", fontWeight = "bold")
  ) #end mdist output
  
  output$mdist_sum2 <- renderDT({summarize_mdist_table(getMdist2()[[1]])})
  
  get_grp2_idx <- reactive({
    nms <- colnames(getMdist2()[[1]])
    return(which(grepl(input$mdist_iden2, nms)))
  })
  
  output$mdist_grp2 <- renderDataTable(datatable(group_mdist_table(getMdist2()[[1]],
                                                                   get_grp2_idx())) %>%
                                         formatStyle("Outcome", fontWeight = "bold") %>%
                                         formatStyle("Outcome", target="row",
                                                     backgroundColor = styleEqual(
                                                       levels = c("Pass","Fail"),
                                                       values = c("#e6ffee","#ffe6e6")))
  ) #end mdist output
  
  #data download handler
  output$data2Download <- downloadHandler(
    filename = "Sample.csv",
    content = function(file) {
      mean_mdist <- getMdist2()[[4]]
      date_reported <- as.character(Sys.Date())
      head <- rbind("Date Reported"=date_reported,"Mean M-Distance"=mean_mdist,"")
      write.table(head, file,sep=",",quote=FALSE,row.names = TRUE,col.names = FALSE)
      write.table(getMdist2()[[1]],file,append=TRUE,sep=",", row.names = FALSE)
      write.table("",file,append=TRUE,sep=",",row.names = FALSE, col.names = FALSE)
      write.table("Quantiles per Outcome",file,append=TRUE,sep=",",
                  row.names = FALSE, col.names = FALSE)
      write.table(as.matrix(t(c("","Passes", "Fails"))),
                  file, append=TRUE,sep=",",row.names = FALSE, col.names = FALSE)
      write.table(summarize_mdist_table(getMdist2()[[1]]),file,append=TRUE,sep=",", row.names = TRUE, col.names = FALSE)
      write.table("",file,append=TRUE,sep=",",row.names = FALSE, col.names = FALSE)
      write.table(paste0("Per ", input$mdist_iden2, " Data"),file,append=TRUE,sep=",", row.names = FALSE, col.names = FALSE)
      write.table(group_mdist_table(getMdist2()[[1]], get_grp2_idx()),file,append=TRUE,sep=",", row.names = FALSE)
    }
    # content = function(file) {
    #   write.csv(getMdist2()[[1]], file, row.names = TRUE)
    # }
  ) #end data download handler
  
  #data download handler (spectral data)
  output$downloadData2 <- downloadHandler(
    filename = function() {
      paste0("sample_spectral_data.csv")
    },
    content = function(file) {
      df <- getSampleData()[[1]]
      scan_sample <- rownames(df)
      df <- cbind(scan_sample, df)
      write.csv(df, file, row.names = FALSE)
    }
  ) #end data download handler
  
  #plot download handler
  output$downloadPlot2 <- downloadHandler(
    filename = function() {
      paste0("sample_spectral_plot.png")
    },
    content = function(file) {
      png(file, width=1000)
      spectraPlot2()
      dev.off()
    }
  ) #end plot download handler
  
  # Sample Normality Plot
  output$sample.normalityplot <- renderPlot(NormalityPlot(getMdist2()[[1]], input$sample.normalityselect,
                                                          tab.name='Sample',kind=input$sample.normalityplottype))
  
  #plot download handler (normality plot)
  output$downloadSampNormPlot <- downloadHandler(
    filename = function() {
      paste0("sample_normality_plot.png")
    },
    content = function(file) {
      ggsave(file, NormalityPlot(getMdist2()[[1]], input$sample.normalityselect,
                                 tab.name='Sample',kind=input$sample.normalityplottype), 
             width=2500, height=1125, units='px')
    }
  ) #end plot download handler
  
  
  #_________________________________________________________________________________________________________________________________
  #GENERIC FUNCTIONS
  
  load_tellspec_file <- function(datapath, fname="") {
    # Open the input filepath

    dataset <- as.data.frame(read.csv2(datapath, header = FALSE, sep = ",", 
                                       fill = TRUE))
    
    # Remove first row on some headerless .csvs
    if(!is.na(dataset[1,1]) & dataset[1,1] == "sep=") {
      dataset <- dataset[-1,]}
    
    # Getting colnames when spectral data titles are ""
    if(any(which(dataset[1,] == ""))) {
      first_dat_col <- which(dataset[1,] == "")[1]
      colnames(dataset) <- c(dataset[1,][1:first_dat_col-1], seq(0,ncol(dataset)-first_dat_col))
      dataset <- dataset[-1,]
    }
    
    # Getting colnames when the first row is entirely filled
    else if (!any(which(dataset[1,] == ""))) {
      colnames(dataset) <- dataset[1,]
      dataset <- dataset[-1,]
    }
    
    # Change column names starting with "..."
    if (any(grepl("...", colnames(dataset), fixed=TRUE))) {
      bns <- length(colnames(dataset)[grepl("...", colnames(dataset), fixed=TRUE)])
      bns <- seq(0,bns)
      colnames(dataset)[grepl("...", colnames(dataset), fixed=TRUE)] <- bns
      
      #dataset <- dataset[-1,]
    }
    
    # Change column names with X/V then a number
    else if (any(grepl("[XV][0-9]*", colnames(dataset)))) {
      if (any(grepl("X[0-9]*", colnames(dataset)))) {pattern <- "X[0-9]*"}
      else {pattern <- "V[0-9]*"}
      
      bns <-length(colnames(dataset)[grepl(pattern, colnames(dataset), fixed=TRUE)])
      bns <- seq(0,bns)
      colnames(dataset)[grepl(pattern, colnames(dataset))] <- bns
    }
    
    # Adds filename before sample name
    fname_addition <- paste0(tools::file_path_sans_ext(fname),"_")
    dataset[,1][seq(2,length(dataset[,1]), 2)] <- paste0(fname_addition, 
                                                         dataset[,1][seq(2,length(dataset[,1]), 2)])
    
    # Change column titles to wavelengths in nm.
    datcols <- grepl("^(?=.)([+-]?([0-9]*)(\\.([0-9]+))?)$", colnames(dataset), perl=TRUE)
    if (dataset[1,1] != "") {colnames(dataset)[datcols] <- dataset[2,datcols]}
    else {colnames(dataset)[datcols] <- dataset[1,datcols]}
    
    return(dataset)
  }
  
  load_tiapp_file <- function(datapath, fname="") {
    # Open the input filepath
    dataset <- as.data.frame(read.csv(datapath,
                                      header = FALSE,
                                      sep = ",",  #Needs to be comma delineated not tab
                                      stringsAsFactors = FALSE))
    
    # Grab Data Elements
    Timestamp <- dataset[2,2]
    System_Temperature <- as.numeric(dataset[4,2]) 
    Detector_Temperature <- as.numeric(dataset[5,2])
    Humidity <- as.numeric(dataset[6,2])
    Sensor_Serial <- dataset[12,2]
    Scan_Config <- dataset[13,2]
    PGA <- as.numeric(dataset[24,2])
    spec_data <- as.numeric(dataset[27:nrow(dataset),4]) / PGA # adjust by PGA
    wv_data <- as.numeric(dataset[27:nrow(dataset),1])
    
    # Wv_row creation
    SampleName <- strsplit(fname,Scan_Config)[[1]][1]
    wv_row <- data.frame(SampleName="",Timestamp="",Sensor_Serial="",
                         Scan_Config="",System_Temperature="",
                         Detector_Temperature="",Humidity="",PGA="", 
                         t(wv_data))
    
    # Data_Row creation
    metadata <- cbind(SampleName, Timestamp,Sensor_Serial,Scan_Config,
                      System_Temperature,Detector_Temperature,Humidity,PGA)
    data_row <- data.frame(metadata, t(spec_data))
    
    # Create two row output dataframe
    Out <- rbind(wv_row, data_row)
    last_meta_colnum <- which(colnames(Out) == "PGA")
    colnames(Out)[(last_meta_colnum+1):ncol(Out)] <- t(wv_data)#seq(0,ncol(Out)-last_meta_colnum-1)
    
    return(Out)
  }
  
  load_compressed_file <- function(datapath, fname = "") {
    # Open the input filepath
    dataset <- as.data.frame(read.csv2(datapath, header = FALSE, sep = ",", 
                                         fill = TRUE))
      
    # Remove first row on some headerless .csvs
    if(!is.na(dataset[1,1]) & dataset[1,1] == "sep=") {
      dataset <- dataset[-1,]
    }
    
    colnames(dataset) <- dataset[1,]
    dataset <- dataset[-1,]
    
    fname_addition <- paste0(tools::file_path_sans_ext(fname),"_")
    dataset[,1] <- paste0(fname_addition, dataset[,1])
    
    return(dataset)
  }
  
  load_labspec_file <- function(datapath, fname = "") {
    df <- as.data.frame(t(read.csv(datapath,
                                   header = FALSE,
                                   sep = "\t",
                                   #skip = 33, 
                                   #row.names = 1, 
                                   stringsAsFactors = FALSE)))
    
    #complete cases
    df <- df[complete.cases(df),]
    
    # Get Data
    wave_colnum <- which(grepl("^(?=.)([+-]?([0-9]*)(\\.([0-9]+))?)$", df[1,], perl=TRUE))[1] - 1
    df <- df[wave_colnum:ncol(df)]
    colnames(df) <- df[1,]
    df <- df[2:nrow(df),]
    
    return(df)
  }
  
  load_all_filetypes <- function(filetype, datapaths, fnames) {
    
    # Tellspec
    if (filetype == "tellspec") {
      dataset <- load_tellspec_file(datapaths[1], fnames[1])
      if (length(datapaths) > 1) {
        for (i in c(2:length(datapaths))) {
          new_file <- load_tellspec_file(datapaths[i], fnames[i])
          dataset <- rbind(dataset, new_file)
        }
      }
    }
    
    # TI App
    else if (filetype == "tiapp") {
      dataset <- load_tiapp_file(datapaths[1], fnames[1])
      if (length(datapaths) > 1) {
        
        # Preallocate space based off number of input files
        preallocated_df <- data.frame(matrix(NA, nrow = 2*length(datapaths), 
                                             ncol = ncol(dataset)))
        preallocated_df[1:2,] <- dataset
        colnames(preallocated_df) <- colnames(dataset)
        j <- 3
        
        # Fill in remaining data
        for (i in c(2:length(datapaths))) {
          new_file <- load_tiapp_file(datapaths[i], fnames[i])
          preallocated_df[j:(j+1),] <- new_file
          j <- j + 2
        }
        
        dataset <- as.data.frame(preallocated_df)
      }
    }
    
    # Compressed
    else if (filetype == "compressed") {
      dataset <- load_compressed_file(datapaths[1], fnames[1])
      
      if (length(datapaths) > 1) {
        for (i in c(2:length(datapaths))) {
          new_file <- load_compressed_file(datapaths[i], fnames[i])
          dataset <- rbind(dataset, new_file)
        }
      }
    }
    
    else if (filetype == "labspec") {
      dataset <- load_labspec_file(datapaths[1], fnames[1])
      
      if (length(datapaths) > 1) {
        for (i in c(2:length(datapaths))) {
          new_file <- load_labspec_file(datapaths[i], fnames[i])
          dataset <- rbind(dataset, new_file)
        }
      }
    }
    
    if (filetype == "tellspec" | filetype == "tiapp") {
      if (dataset[1,1] != "") {
        dataset <- dataset[seq(1,nrow(dataset),2),]
      }
      else {dataset <- dataset[seq(2,nrow(dataset),2),]}
    }
    
    dataset <- dataset %>% mutate_at(colnames(dataset)[grepl("^(?=.)([+-]?([0-9]*)(\\.([0-9]+))?)$", 
                                                             colnames(dataset), perl=TRUE)], as.numeric)
    
    return(dataset)
  }
  
  group_mdist_table <- function(mdist_table, avg_by_idx) {
    
    # Collect colnames to specify correct operations
    num_cols <- c("M-Dist","Residual","Factors","Cutoff","Difference")
    mdist_table <- mdist_table %>% mutate_at(num_cols, as.numeric)
    str_cols <- c(colnames(mdist_table)[!(colnames(mdist_table) %in% num_cols)])
    
    # Split mdist_table and add column to group relative to
    Col_To_Avg_By <- mdist_table[avg_by_idx]
    colnames(Col_To_Avg_By) <- "Avgby"
    mdist_table$Avgby <- Col_To_Avg_By
    df_split1 <- mdist_table[c("Avgby", num_cols)]
    df_split2 <- mdist_table[c("Avgby", str_cols)]

    # Average columns with numerical data
    df_split1 <- df_split1 %>% group_by(Avgby) %>% summarise(across(all_of(num_cols), mean))
    
    # Summarize Outcome column and collapse remaining columns relative to column to group
    outcme <- df_split2 %>% group_by(Avgby) %>% summarise(across("Outcome", n_distinct))
    df_split2 <- df_split2 %>% group_by(Avgby) %>% summarise(across(all_of(str_cols), first))
    df_split2[outcme$Outcome > 1,"Outcome"] <- "Disagreements Found"
    df_split2[outcme$Outcome == 0,"Outcome"] <- "No Data Found"
    
    # Merge split dataframes
    mdist_table <- cbind(df_split2[2], df_split1[2:ncol(df_split1)], df_split2[3:ncol(df_split2)])
    
    return(mdist_table)
  }
  
  summarize_mdist_table <- function(mdist_table) {
    passes <- as.numeric(mdist_table[mdist_table$Outcome == "Pass",][,2])
    fails <- as.numeric(mdist_table[mdist_table$Outcome == "Fail",][,2])
    
    num_samples <- c(length(passes), length(fails))
    pass_quarts <- quantile(passes, probs = c(0,0.25,0.5,0.75,1))
    fail_quarts <- quantile(fails, probs = c(0,0.25,0.5,0.75,1))
    
    sum_table <- cbind(pass_quarts, fail_quarts)
    sum_table <- rbind(num_samples, sum_table)
    
    rownames(sum_table) <- c("Number of Samples", "Minimum", "25% Quartile",
                             "Mean", "75% Quartile", "Maximum")
    colnames(sum_table) <- c("Passes", "Fails")
    
    sum_table <- as.data.frame(sum_table)
    return(sum_table)
  }
  
  truncate_by_wavelength <- function(spectra, start_nm, end_nm) {
    lowerbound_idx <- which(as.numeric(colnames(spectra)) >= start_nm)[1]
    upperbound_idx <- which(as.numeric(colnames(spectra)) <= end_nm)[length(which(as.numeric(colnames(spectra)) <= end_nm))]
    
    return(spectra[lowerbound_idx:upperbound_idx])
  }
  
  auto_truncate <- function(ref_data, samp_data) {
    # Get colnames (wavelengths) for ref and samp data
    colnames_ref <- colnames(ref_data)
    colnames_samp <- colnames(samp_data)
    
    # Get the first and final reference wavelengths
    min_wavelength_ref <- colnames_ref[1]
    max_wavelength_ref <- colnames_ref[length(colnames_ref)]
    
    # Get First wavelength index
    if (colnames_ref[1] %in% colnames_samp) {
      samp_start <- which(colnames_samp == colnames_ref[1])
    }
    else {
      ref_floor <- floor(as.numeric(colnames_ref[1]))
      samp_start <- which(as.numeric(colnames_samp) > ref_floor)[1]
    }
    
    # Get final wavelength index
    if (colnames_ref[length(colnames_ref)] %in% colnames_samp) {
      samp_end <- which(colnames_samp == colnames_ref[length(colnames_ref)])
    }
    else {
      ref_ceil <- ceiling(as.numeric(colnames_ref[length(colnames_ref)]))
      samp_end <- which(as.numeric(colnames_samp) < ref_ceil)[length(which(as.numeric(colnames_samp) < ref_ceil))]
    }
    
    # Adjust wavelength indexes to fit number of ref columns
    ## sample truncation too small
    while (length(colnames_ref) > length(colnames_samp[samp_start:samp_end])) {
      start_diff <- abs(as.numeric(colnames_samp[samp_start]) - as.numeric(colnames_ref[1]))
      end_diff <- abs(as.numeric(colnames_samp[samp_end]) - as.numeric(colnames_ref[length(colnames_ref)]))
      
      if (start_diff > end_diff) { # Expand sample trunc towards most difference from reference wavelength
        if (samp_start > 1) {
          samp_start <- samp_start - 1 # Decrease first index, unless we reach the start
        }
        else {break}
      }
      
      else {
        if (samp_end < length(colnames_samp)) {
          samp_end <- samp_end + 1 # Increase final index, unless we reach the end
        }
        else {break}
      }
    }
    
    ## sample truncation too large
    while (length(colnames_ref) < length(colnames_samp[samp_start:samp_end])) {
      start_diff <- abs(as.numeric(colnames_samp[samp_start]) - as.numeric(colnames_ref[1]))
      end_diff <- abs(as.numeric(colnames_samp[samp_end]) - as.numeric(colnames_ref[length(colnames_ref)]))
      
      if (start_diff > end_diff) { # Shrink sample trunc towards most difference from reference wavelength
        if (samp_start < samp_end) {
          samp_start <- samp_start + 1 # Increase first index, unless we reach start index
        }
        else {break}
      }
      
      else {
        if (samp_end > samp_start) {
          samp_end <- samp_end - 1 # Decrease final index, unless we reach final index
        }
        else {break}
      }
    }
    
    return(samp_data[samp_start:samp_end])
  }
  
  NormalityPlot <- function(df, target.col, tab.name='Reference', kind='PDF') {
    colnames(df)[which(colnames(df) == 'M-Dist')] <- "M.Dist"
    df <- df %>% mutate_at("M.Dist", as.numeric)
    df["Target"] <- df[,target.col]
    plt.title <- paste(tab.name, "M-Distance Distributions by", target.col)
    
    plt <- ggplot(data=df, aes(x=M.Dist, y=Target, fill=Target))
    
    if(kind=='PDF') { # PDF
      plt <- plt + stat_halfeye(fill_type='segments', normalize='groups', alpha=0.5,
                                .width=c(0.5, 0.95), trim=FALSE, expand=TRUE)
    }
    else if(kind=='CDF') { # CDF
      plt <- plt + stat_cdfinterval(fill_type='segments', normalize='groups', alpha=0.5,
                                .width=c(0.5, 0.95), trim=TRUE, expand=FALSE,
                                density='bounded', outline_bars=TRUE)
        
      #plt <- ggplot(data=df, aes(x=stage(M.Dist, after_stat=cdf), y=Target, fill=Target))
    }
    
    plt <- plt + stat_summary(geom="point", fun=median) +
      geom_vline(xintercept = getMdist()[[9]], col = "grey30", lty = "dashed") +
      labs(x="M-Distance",y=target.col, fill=target.col, title=plt.title) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme_bw()
    
    return(plt)
  }
}
