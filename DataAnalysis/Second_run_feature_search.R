# parallel feature number search
# U19 pilot data analysis
# Date: 0902022
# created by Ray Su
# Data source: Dr. Qiu Chuan
# Data preparation scrip: DataPrep_BoneStudy119.R
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library("tidyverse")
library("dplyr")
library("R.utils") # Utility functions useful when programming and developing R packages.
library("mvtnorm")
library("tester") # tester allows you to test characteristics of common R objects.
library("turner") # Package designed for working with vectors and lists of vectors, mainly for turning them into other indexed data structures.
# library("diagram") # Functions for Visualizing Simple Graphs (Networks), Plotting Flow Diagrams
# library("shape") # Functions for plotting graphical shapes such as ellipses, circles, cylinders, arrows, ...
library("amap")  # Tools for Clustering and Principal Component Analysis (With robust methods, and paralleled functions).
library(nnet) # multinomial regression
require("parallel")
require("doParallel")
require("foreach")
rm(list = ls())

load("Rcode/RealDataAnalysis/LOS_U19_pilot.RData")
rm(list=c("GeneExpression","RRBS","metabo","tmp.data"))
run1 <- readRDS(file = "Rcode/RealDataAnalysis/Final_LOS_results.RDS")
feature_slection <- run1$model$outer_model %>% 
  filter(weight!=0) %>% 
  transform(name=as.character(name)) %>% 
  filter(name!="BMD")
# View(head(dfCompile))
dfCompile <- dfCompile %>% select(-feature_slection$name)
# Import libraries
# library("tidyverse")

# doParallel::registerDoParallel(cl)
cl <- makeCluster(40)
doParallel::registerDoParallel(cl)





sourceDirectory("2018_msPLS-master/Example/sPLSPM/R/")
files.sources = list.files(path = "2018_msPLS-master/Example/sPLSPM/R/",full.names = T)
sapply(files.sources, source)

# colnames(dfCompile)[1:3]
# View(dfCompile[,1:5])
# Gene expression X1
EXPL_X1 = c(0,1,0,0)
# Methylation X2
EXPL_X2 = c(0,0,1,0)
# metabolites X3
EXPL_X3 = c(1,1,0,0)
# BMD (X4)
RESP_Y = c(1,1,1,0)
path_matrix = rbind(EXPL_X1, EXPL_X2,EXPL_X3,RESP_Y)

# blocks of outer model
blocks <- list(grep(pattern = "GEID*",x = colnames(dfCompile),value = F),
               grep(pattern = "RRBSID*",x = colnames(dfCompile),value = F),
               grep(pattern = "metaID*",x = colnames(dfCompile),value = F),
               grep(pattern = "BMD",x = colnames(dfCompile)))
modes = c("B","B","B","A")

dfCompile <- as.data.frame(dfCompile)
# lapply(blocks,function(x) length(x)*0.01)

tmp.select <- seq(from=30,to=150,by=10)
gd.table <- expand.grid(tmp.select,tmp.select,tmp.select)
# iter=1
# as.vector(t(gd.table[iter,]))

myRDS <- list.files(path = "Rcode/RealDataAnalysis/FeatureNumberSearch/run2",pattern = ".RDS")
complet.run <- word(string = myRDS,start = 3,end = 3,sep = "_") %>% gsub(pattern = ".RDS",replacement = "")
continu.run <- c(1:nrow(gd.table))[-as.numeric(complet.run)]
cat("statting jobs \n")
# iter=1682
# foreach(iter=1:nrow(gd.table),.packages = c("tester","turner","mvtnorm","amap","nnet")) %dopar% {
# for (iter in continu.run) {
foreach(iter=continu.run,.packages = c("tester","turner","mvtnorm","amap","nnet")) %dopar% {
  ## fit the model
  time_data <- system.time(
    s_satpls <- splspm(dfCompile, path_matrix, blocks, modes, scheme="path",
                       scaled=F, penalization = "enet", nonzero = as.vector(t(gd.table[iter,])),
                       # scaled=F, penalization = "enet",
                       cross_validate = F, lambda = 1, maxiter = 500)
  )
  
  ############### msPLSDA ##############
  # X1 = gene expression
  # X2 = Methylome
  # X3 = Metabolites
  
  my_score_data <- data.frame(X1=s_satpls$scores[,1],
                              X2=s_satpls$scores[,2],
                              X3=s_satpls$scores[,3])
  
  my_score_data <- cbind(my_score_data,dfCompile$BMD)
  colnames(my_score_data)[4] <- "myGroup"
  
  # my_score_data$myGroup <- relevel(my_score_data$myGroup, ref = "1")
  training_model <- multinom(myGroup ~ X1 + X2 + X3, data = my_score_data)
  training_outcom_truth <- as.factor(my_score_data$myGroup)
  pp <- fitted(training_model)
  # Get the label based on the fitted prob.
  # max_prob_position <- apply(pp, 1, which.is.max)
  max_prob_position <- ifelse(pp>0.5,1,0)
  max_prob_position <- factor(max_prob_position)
  performance.output.training.overall <- caret::confusionMatrix(max_prob_position,training_outcom_truth)
  tmp.filename <- paste("Rcode/RealDataAnalysis/FeatureNumberSearch/run2/Out_file_",iter,".RDS",sep = "")
  mylist <- list("Selection"=as.vector(t(gd.table[iter,])),
                 "PredictionOut"=performance.output.training.overall,
                 "model"=s_satpls)
  # saveRDS(object = mylist,file = tmp.filename)
}

# stopCluster()
on.exit(stopCluster(cl))


