# Run MOFA feature selection in simulation study
# https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/downstream_analysis.html
# Created by Ray Su
# Date: 12282021
# This script is used to perform MOFA on the simulated data 
# msPLSDA project
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("MOFA2")

library("tidyverse")
require("parallel")
require("doParallel")
# library("doMC")
require("foreach")
library("MOFA2")
library("data.table")

# registerDoParallel(25)
# MOFA method
# Load simulated Data 
rm(list=ls())
# registerDoParallel(30)

# all.filepath <- list.dirs(path = "SimData")[-1]
all.filepath <- paste("SimData/SimData_",750,sep = "")
for (mydir in all.filepath) {
  # mydir <- list.dirs(path = "SimData")[5]
  
  print(mydir)
  mydata <- list.files(path = mydir , pattern = "*.RDS",full.names = T)
  runID <- word(string =mydir ,2,sep = "_")
  print(runID)
  
  my.bind.out <- list()
  # my.bind.out <- foreach(iter=1:length(mydata),.packages = c("tidyr","dplyr","MOFA2")) %dopar% {
  # my.bind.out <- foreach(iter=1:3,.packages = c("tidyr","dplyr","MOFA2"),.errorhandling = 'pass') %dopar% {
  for( iter  in 1:length(mydata)){
  # for( iter  in 1:3){
    # iter=2
    tmp.iter <- mydata[iter]
    my.seed <- word(string = tmp.iter,4,sep = "_") %>% gsub(pattern = ".RDS",replacement = "",x = .)
    # print(my.seed)
    set.seed(my.seed)
    tmp.input <- readRDS(file = tmp.iter)
    Data <- tmp.input[[1]] 
    
    parameters <- tmp.input[[2]]
    
    # blocks for the omics datasets
    blocks <- list("data1"=(grep(pattern = "*or1_*",x = colnames(Data),value = F)),
                   "data2"=(grep(pattern = "*or2_*",x = colnames(Data),value = F)))
    
    data_block <- list("data1"=t(Data[,blocks[[1]]]),
                        "data2"=t(Data[,blocks[[2]]]))
    
    colnames(data_block[[1]]) <- paste("sample_",1:nrow(Data),sep="")
    colnames(data_block[[2]]) <- paste("sample_",1:nrow(Data),sep="")
    
    MOFAobject <- create_mofa(data = data_block,groups = Data$GroupInfo)
    # MOFAobject
    # plot_data_overview(MOFAobject)
    # Define data options
    data_opts <- get_default_data_options(MOFAobject)
    # Define model options
    model_opts <- get_default_model_options(MOFAobject)
    model_opts$num_factors <- 2
    
    train_opts <- get_default_training_options(MOFAobject)
    # train_opts$save_interrupted <- TRUE
    # train_opts$outfile <- paste("Results/MOFA/MOFA_",runID,"_",my.seed,".hdf5",sep = "")
    train_opts$seed <- my.seed
    
    MOFAobject <- prepare_mofa(
      object = MOFAobject,
      data_options = data_opts,
      model_options = model_opts,
      training_options = train_opts
    )
    
    # MOFAobject.trained <- run_mofa(MOFAobject, outfile)
    MOFAobject.trained <- run_mofa(MOFAobject,
                                   # outfile =  paste("Results/MOFA/MOFA_",runID,"_",my.seed,".hdf5",sep = ""),
                                   save_data = F,use_basilisk = TRUE)
    # plot_variance_explained(MOFAobject.trained, x="view", y="factor")
    # plot_weights(MOFAobject.trained,
    #              view = "data1",
    #              factor = 1,
    #              nfeatures = 50,     # Number of features to highlight
    #              scale = T,          # Scale weights from -1 to 1
    #              abs = F             # Take the absolute value?
    # )
    # plot_top_weights(MOFAobject.trained,
    #                  view = "data1",
    #                  factor = 1,
    #                  nfeatures = parameters$SelectedFeatures
    # )
    # cat("run:",runID,"go go!\n")
    selection.features <- get_weights(MOFAobject.trained,scale = T,abs = T,as.data.frame = T) %>% 
      filter(factor=="Factor1") %>% 
      group_by(view) %>% 
      slice_max(value,n = parameters$SelectedFeatures)
    
    selection.d1 <- selection.features %>% 
      filter(view=="data1") %>% 
      dplyr::select(feature) %>% pull %>% 
      gsub(pattern = "\\d_\\d+",replacement = "",x = .) %>% 
      table()
    
    selection.d2 <- selection.features %>% 
      filter(view=="data2") %>% 
      dplyr::select(feature) %>% pull %>% 
      gsub(pattern = "\\d_\\d+",replacement = "",x = .) %>% 
      table()
    
    names(selection.d1) <- paste(names(selection.d1),"_EXPL_X",sep = "")
    names(selection.d2) <- paste(names(selection.d2),"_EXPL_Y",sep = "")
    feature.count <- c(selection.d1,selection.d2)
    # Summarize variable selections for each block
    
    # feature.count <- feature.count[-c(1,2)]
    nms <-  c("discrCor_EXPL_X", "discrUncor_EXPL_X", "noiseCor_EXPL_X",
              "noiseUncor_EXPL_X", "discrCor_EXPL_Y", "discrUncor_EXPL_Y",
              "noiseCor_EXPL_Y", "noiseUncor_EXPL_Y")
    Missing <- setdiff(nms,names(feature.count))
    feature.count[Missing] <- 0
    # Get the results of feature selection count 
    feature.count <- feature.count[nms]
    
    # Gather all results
    output.return <- list("feature.count" =feature.count,
                          "seed"=my.seed,
                          "parameters"=parameters)
    my.bind.out[[iter]] <- output.return
    # return(output.return)
  }
  # for(iter in 1:length(mydata)){
  # for(iter in 1:3){
  # print(iter)
  # iter <- 250
  
  saveRDS(object = my.bind.out,file = paste("Results/MOFA/MOFA_",runID,".RDS",sep = ""))
  cat("Output Saved:", paste("Results/MOFA/MOFA_",runID,".RDS",sep = ""), "\n")
}

# stopImplicitCluster()

