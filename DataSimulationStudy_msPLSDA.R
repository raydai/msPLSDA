#Workable
rm(list=ls())
# library("tidyverse")
# library("sPLSPM")
library("tidyr")
library("dplyr")
library("holodeck")
# library("fastDummies")   # convert the groups into dummy variables
#parameters for testing
# library(sPLSPM)
library(R.utils) # Utility functions useful when programming and developing R packages.
library(mvtnorm)
library(tester) # tester allows you to test characteristics of common R objects.
library(turner) # Package designed for working with vectors and lists of vectors, mainly for turning them into other indexed data structures.
library(diagram) # Functions for Visualizing Simple Graphs (Networks), Plotting Flow Diagrams
library(shape) # Functions for plotting graphical shapes such as ellipses, circles, cylinders, arrows, ...
library(amap)  # Tools for Clustering and Principal Component Analysis (With robust methods, and paralleled functions).
require(parallel)
require(doParallel)
require(foreach)
setwd("/lustre/project/hdeng2/ksu2/msPLSDA/msPLSDA/EnetPLS/")

sourceDirectory("/lustre/project/hdeng2/ksu2/msPLSDA/msPLSDA/EnetPLS/2018_msPLS-master/Example/sPLSPM/R")
files.sources = list.files(path = "/lustre/project/hdeng2/ksu2/msPLSDA/msPLSDA/EnetPLS/2018_msPLS-master/Example/sPLSPM/R/",full.names = T)
sapply(files.sources, source)

# args <- c(2,1,2,3,50,0.01,50,200,30,168,1,1)
args<-commandArgs(TRUE)
# follow the para-table order

# ID
ID <- as.integer(args[2])
# OmicSets
OmicSets <- as.integer(args[3])
# CateOutcome
CateOutcome <- as.integer(args[4])
# SampleSize
SampleSize <- as.integer(args[5])
# Sparsity
Sparsity <- as.numeric(args[6])
# SelectedFeatures
SelectedFeatures <- as.integer(args[7])
# TotalFeatureInOmics
TotalFeatureInOmics <- as.numeric(args[8])
# noiseCor
noiseCor <- as.numeric(args[9])
# noiseUncor
noiseUncor <- as.numeric(args[10])
# discrCor
discrCor <- as.numeric(args[11])
# discrUncor
discrUncor <- as.numeric(args[12])
# batch id
tmp.round <- as.numeric(args[13])



tmp.noisy.cor.proportion <- 0.15
tmp.discr.cor.proportion <- 0.5
tmp.mean = c(-0.3,-0.2,0.3) # the mean for three groups
#### Testing variables
# Create three groups
# SampleSize=80
# OmicSets=2
# SelectedFeatures =50
# CateOutcome= 3 # The number of categorical outcomes

# 
# # all features for each modelity
# TotalFeatureInOmics=500
# # we may call this dense instead of sparsity! 
# Sparsity = 0.1

# # noise data
# tmp.num.noise <- ceiling(TotalFeatureInOmics*(1-Sparsity))
# noiseCor <- ceiling(tmp.num.noise*tmp.noisy.cor.proportion)
# noiseUncor <-  tmp.num.noise-noiseCor
# # discriminators features

# tmp.discr <- TotalFeatureInOmics - tmp.num.noise
# discrCor <- ceiling(tmp.discr*tmp.discr.cor.proportion)
# discrUncor <- tmp.discr-discrCor

# msPLSDA_performance(SampleSize = 80,CateOutcome = 3,discrCor = discrCor,discrUncor = discrUncor,noiseCor = noiseCor,noiseUncor = noiseUncor,tmp.mean = tmp.mean,seed.num = 1,SelectedFeatures =SelectedFeatures,OmicSets =  OmicSets)

msPLSDA_performance <- function(SampleSize,CateOutcome,discrCor,discrUncor,noiseCor,noiseUncor,tmp.mean,seed.num,SelectedFeatures,OmicSets) {
  set.seed(seed.num)
  # Generate data
  tmp.data <- sim_cat(n_obs = SampleSize,n_groups = CateOutcome,name = "GroupInfo") %>% 
    group_by(GroupInfo) %>% 
    sim_discr(n_vars = discrCor, var = 1,cov = 0.7,group_means = tmp.mean,name = "discrCor1") %>% 
    sim_discr(n_vars = discrUncor, var = 1,cov = 0.0,group_means = tmp.mean,name = "discrUncor1") %>%
    ungroup() %>% 
    sim_covar(n_vars = noiseCor, var = 1, cov = 0.7, name = "noiseCor1") %>% 
    sim_covar(n_vars = noiseUncor, var = 1, cov = 0.0, name = "noiseUncor1") %>% 
    group_by(GroupInfo) %>% 
    sim_discr(n_vars = discrCor, var = 1,cov = 0.3,group_means = tmp.mean,name = "discrCor2") %>% 
    sim_discr(n_vars = discrUncor, var = 1,cov = 0.0,group_means = tmp.mean,name = "discrUncor2") %>%
    ungroup() %>% 
    sim_covar(n_vars = noiseCor, var = 1, cov = 0.3, name = "noiseCor2") %>% 
    sim_covar(n_vars = noiseUncor, var = 1, cov = 0.0, name = "noiseUncor2") 
  
  # Create a temporal for random sampling
  tmp.data <- tmp.data%>% 
    mutate(ID=seq_len(nrow(tmp.data)))
  
  # transform(myGroup=as.factor(myGroup))
  # plot heatmapping for the simulation datasets
  # tmp.data %>%select(-c("GroupInfo","ID")) %>% 
  #   cov() %>% 
  #   heatmap(Rowv = NA,Colv = NA,symm = T)
  
  # library(dplyr)
  training_IDs <- tmp.data %>% group_by(GroupInfo) %>% sample_n(ceiling(SampleSize*0.7/CateOutcome)) %>% select("ID") %>% pull
  # length(training_IDs)
  # split the tmp.data into training (70%) and testing (70%) datasets 
  '%!in%' <- function(x,y) !(x%in%y)
  # training_df2 <- tmp.data %>% filter(ID%in%training_IDs) %>% select(-ID) %>%
  #   transform(GroupInfo=as.factor(GroupInfo))
  training_df <- tmp.data %>% filter(ID%in%training_IDs) %>% select(-ID) %>%
    transform(GroupInfo=ifelse(GroupInfo=="a",1,ifelse(GroupInfo=="b",2,3)))
  
  # training_df <- tmp.data %>% filter(ID%in%training_IDs) %>% select(-ID) %>%
  #   fastDummies::dummy_cols(., select_columns = "GroupInfo") %>%
  #   select(-GroupInfo)
  # View(training_df[,c(1001:1003)])
  # training_df <- tmp.data %>% filter(ID%in%training_IDs) %>% select(-ID) %>%%>%
  #   transform(GroupInfo=ifelse(GroupInfo=="a",1,ifelse(GroupInfo=="b",2,3)))
  # 
  testing_df <- tmp.data %>% filter(ID%!in%training_IDs) %>% select(-ID)%>%
    transform(GroupInfo=ifelse(GroupInfo=="a",1,ifelse(GroupInfo=="b",2,3)))
  
  
  # sourceDirectory("2018_msPLS-master/Example/sPLSPM/R/")
  # files.sources = list.files(path = "2018_msPLS-master/Example/sPLSPM/R/",full.names = T)
  # sapply(files.sources, source)
  
  
  # Assign training dataset
  Data <- training_df
  # Data <- cbind(X,Y,Z)
  EXPL_X = c(0,1,0)
  EXPL_Y = c(1,0,0)
  RESP_Z = c(1,1,0)
  path_matrix = rbind(EXPL_X, EXPL_Y,RESP_Z)
  
  # blocks of outer model
  blocks <- list((grep(pattern = "*or1_*",x = colnames(Data),value = F)),
                 (grep(pattern = "*or2_*",x = colnames(Data),value = F)),
                 grep(pattern = "GroupInfo",x = colnames(Data)))
  modes = c("B","B","A")
  
  ## fit the model
  time_data <- system.time(
    s_satpls <- splspm(Data, path_matrix, blocks, modes, scheme="path",
                       scaled=F, penalization = "enet", nonzero = rep(x = SelectedFeatures,OmicSets),
                       cross_validate = T, lambda = 1, maxiter = 100)
  )
  # time_data <- system.time(
  #   s_satpls <- splspm(Data, path_matrix, blocks, modes, scheme="path",
  #                      scaled=F, penalization = "enet", nonzero = c(150,150),
  #                      cross_validate = F, lambda = 1, maxiter = 100)
  # )
  # time_data <- system.time(
  #   s_satpls <- splspm(Data, path_matrix, blocks, modes, scheme="path",
  #                      scaled=F, penalization = "ust",
  #                      cross_validate = F, lambda = 1, maxiter = 500)
  # )
  # summary(s_satpls)
  
  
  ########################################################################################################
  
  ########## feature selection ##########
  # Summarize variable selections for each block
  feature_slection <- s_satpls$outer_model %>% 
    filter(weight!=0) %>% 
    mutate("predCate"=gsub(pattern = "\\d_\\d+",replacement = "",x = name)) %>% 
    filter(predCate%in%c("noiseCor","noiseUncor","discrUncor","discrCor")) %>% 
    group_by(block,predCate) %>% summarise(Count=n()) %>% 
    mutate(predCate= paste(predCate,block,sep = "_"),.keep="unused")
  feature_slection <- feature_slection[,-1]
  feature.count <- feature_slection$Count
  names(feature.count) <- feature_slection$predCate
  # feature.count <- feature.count[-c(1,2)]
  nms <-  c("discrCor_EXPL_X", "discrUncor_EXPL_X", "noiseCor_EXPL_X",
            "noiseUncor_EXPL_X", "discrCor_EXPL_Y", "discrUncor_EXPL_Y",
            "noiseCor_EXPL_Y", "noiseUncor_EXPL_Y")
  Missing <- setdiff(nms,names(feature.count))
  feature.count[Missing] <- 0
  # Get the results of feature selection count 
  feature.count <- feature.count[nms]
  
  ############### msPLSDA ##############
  # save.image(file = "test_score.RData")
  # library(cluster)    # clustering algorithms
  # library(factoextra) # clustering algorithms & visualization
  my_score_data <- data.frame(X1=s_satpls$scores[,1],
                              X2=s_satpls$scores[,2])
  
  my_score_data <- cbind(my_score_data,Data$GroupInfo)
  colnames(my_score_data)[3] <- "myGroup"
  
  # km.res <- kmeans(x = my_score_data,centers = 3)
  # fviz_cluster(km.res, data = my_score_data[,c(1,2)])
  library(nnet) # multinomial regression
  # my_score_data$myGroup <- relevel(my_score_data$myGroup, ref = "1")
  training_model <- multinom(myGroup ~ X1 + X2, data = my_score_data)
  training_outcom_truth <- as.factor(my_score_data$myGroup)
  # summary(training_model)
  pp <- fitted(training_model)
  # Get the label based on the fitted prob.
  max_prob_position <- apply(pp, 1, which.is.max)
  max_prob_position <- factor(max_prob_position)
  performance.output.training.overall <- caret::confusionMatrix(max_prob_position,training_outcom_truth)
  
  # max_prob_label <- dplyr::recode(max_prob_position, "1"="A" , "2" = "B" , "3" ="C")
  # table(max_prob_label)
  # sum(table(max_prob_label))
  
  # Calculate the scores for testing dataset
  # s_satpls$outer_model$block %>% unique()
  # get the prediction probability
  # Get X1 score and X2 score
  fit.outer.weight <- s_satpls$outer_model %>% 
    filter(block%in%c("EXPL_X","EXPL_Y")) %>% 
    select(name,weight)
  
  tmp.testing.df <- testing_df %>% select(-GroupInfo) %>% t() %>% as.data.frame()
  tmp.testing.df$name <- rownames(tmp.testing.df)
  tmp.df <- merge(x = tmp.testing.df,y = fit.outer.weight,by = "name")
  rownames(tmp.df) <- tmp.df$name
  tmp.df <- tmp.df %>% select(-name)
  tmp.weight <- tmp.df$weight
  tmp.df <- tmp.df %>% select(-weight) %>% as.matrix()
  tmp.df <- tmp.df*tmp.weight
  # matrix(data = c(1:20),nrow = 2)*c(2,5)
  x1.var <- grep(pattern = "1_",x = rownames(tmp.df))
  tmp.X1 <- tmp.df[x1.var,] %>% colSums()
  tmp.X2 <- tmp.df[-x1.var,] %>% colSums()
  tmp.testing.scores <- data.frame("X1"=tmp.X1,"X2"=tmp.X2)
  testing.predict <- predict(training_model, newdata = tmp.testing.scores, "probs")
  testing.predict <- apply(testing.predict, 1, which.is.max)
  xtab <- table(testing.predict,testing_df$GroupInfo)
  
  fac.level <- unique(testing_df$GroupInfo) %>% sort()
  testing.predict <- factor(x = testing.predict,levels = fac.level)
  testing_truth <- factor(x = testing_df$GroupInfo,levels = fac.level)
  performance.output.testing <- caret::confusionMatrix(testing.predict,testing_truth)
  performance.output.testing.overall <- performance.output.testing$overall
  
  # Gather all results
  output.return <- list("feature.count" =feature.count,
                        "confunsion.testing"=performance.output.testing.overall,
                        "confunsion.training"=performance.output.training.overall,
                        "seed"=my.seed)
  return(output.return)
}

set.seed(446646)
# Goal: 500
iter= 500
my.seed <-sample(x = c(1:10000000),size = iter,replace = F)
# tmp.round <- 2
my.round.seed <- switch(tmp.round,
                        my.seed[1:50],
                        my.seed[50:100],
                        my.seed[101:150],
                        my.seed[151:200],
                        my.seed[201:250],
                        my.seed[251:300],
                        my.seed[301:350],
                        my.seed[351:400],
                        my.seed[401:450],
                        my.seed[451:500])
# The number of thread

# registerDoParallel(cores=(as.integer(args[1])))
# cluster <- makeCluster(detectCores())

NumCores <- detectCores()
cat("The number of cores for mapping: ",NumCores,"\n")
registerDoParallel(NumCores)

# cluster <- makeCluster(as.integer(args[1]))
# registerDoParallel(cluster)
start <- proc.time()
parallel.out <- foreach(i=my.round.seed,.packages = c("tidyr","dplyr","holodeck","tester","turner")) %dopar% {
  msPLSDA_performance(SampleSize = SampleSize,
                      CateOutcome = CateOutcome,
                      discrCor = discrCor,
                      discrUncor = discrUncor,
                      noiseCor = noiseCor,
                      noiseUncor = noiseUncor,
                      tmp.mean = tmp.mean,
                      seed.num = i,
                      SelectedFeatures = SelectedFeatures,
                      OmicSets = OmicSets)  
}
dopar_loop <- proc.time()-start
cat("Running time:", dopar_loop, "\n")
# parallel.out
# stopCluster(cluster)
# Clean up the cluster
stopImplicitCluster()

all.parameters <- list("sample"=SampleSize,
                       "groups" = CateOutcome,
                       "OmicSet"=OmicSets,
                       "Features"=c(discrCor,discrUncor,noiseCor,noiseUncor),
                       "Group.means"= tmp.mean,
                       "Sparsity"=Sparsity,
                       "Desc.proportion"=c(tmp.noisy.cor.proportion,tmp.discr.cor.proportion))    
names(all.parameters$Features) <- c("discrCor","discrUncor","noiseCor","noiseUncor")
names(all.parameters$Desc.proportion) <- c("noise","descriminator")
# Export the simulation results
all.output <- list("Simulation" = parallel.out , "parameter" =all.parameters)
# all.output
dir.create(file.path("./tmp.out_V3"), showWarnings = FALSE)
file.name = paste("./tmp.out_V3/V3_",ID,"_",tmp.round,".rds",sep = "")
saveRDS(all.output, file = file.name)
cat("####Simulation Done####")
# use readRDS(file=file.name) to import output next time.
# For example: 
# mydata <- readRDS(file = "tmp.out/Test.rds")

# mean(as.data.frame(parallel.out)[,1])
