# U19 pilot data analysis
# Date: 08292022
# created by Ray Su
# Data source: Dr. Qiu Chuan
# Data preparation scrip: DataPrep_BoneStudy119.R
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

load("Rcode/RealDataAnalysis/LOS_U19_pilot.RData")
rm(list=c("GeneExpression","RRBS","metabo","tmp.data"))

# Import libraries
# library("tidyverse")
library("dplyr")
library("R.utils") # Utility functions useful when programming and developing R packages.
library("mvtnorm")
library("tester") # tester allows you to test characteristics of common R objects.
library("turner") # Package designed for working with vectors and lists of vectors, mainly for turning them into other indexed data structures.
# library("diagram") # Functions for Visualizing Simple Graphs (Networks), Plotting Flow Diagrams
# library("shape") # Functions for plotting graphical shapes such as ellipses, circles, cylinders, arrows, ...
library("amap")  # Tools for Clustering and Principal Component Analysis (With robust methods, and paralleled functions).
require("parallel")
require("doParallel")
require("foreach")

# setwd("/lustre/project/hdeng2/ksu2/msPLSDA/msPLSDA/EnetPLS/")
# sourceDirectory("/lustre/project/hdeng2/ksu2/msPLSDA/msPLSDA/EnetPLS/2018_msPLS-master/Example/sPLSPM/R")
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

# get the top 10% of features
dfCompile <- as.data.frame(dfCompile)
# lapply(blocks,function(x) length(x)*0.01)
SelectedFeatures= 50
OmicSets= 3
rep(x = SelectedFeatures,OmicSets)

tmp.select <- seq(from=30,to=200,by=10)
gd.table <- expand.grid(tmp.select,tmp.select,tmp.select)

## fit the model
time_data <- system.time(
  s_satpls <- splspm(dfCompile, path_matrix, blocks, modes, scheme="path",
                     scaled=F, penalization = "enet", nonzero = rep(x = SelectedFeatures,OmicSets),
                     # scaled=F, penalization = "enet",
                     cross_validate = F, lambda = 1, maxiter = 500)
)

#An object of class "plspm".
# 
# outer_model	
# Results of the outer model. Includes: outer weights, standardized loadings, communalities, and redundancies
# inner_model	
# Results of the inner (structural) model. Includes: path coeffs and R-squared for each endogenous latent variable
# 
# scores	
# Matrix of latent variables used to estimate the inner model. If scaled=FALSE then scores are latent variables calculated with the original data (non-stardardized).
# 
# path_coefs	
# Matrix of path coefficients (this matrix has a similar form as path_matrix)
# 
# crossloadings	
# Correlations between the latent variables and the manifest variables (also called crossloadings)
# 
# inner_summary	
# Summarized results of the inner model. Includes: type of LV, type of measurement, number of indicators, R-squared, average communality, average redundancy, and average variance extracted
# 
# effects	
# Path effects of the structural relationships. Includes: direct, indirect, and total effects
# 
# unidim	
# Results for checking the unidimensionality of blocks (These results are only meaningful for reflective blocks)
# 
# gof	
# Goodness-of-Fit index
# 
# data	
# Data matrix containing the manifest variables used in the model. Only available when dataset=TRUE
# 
# boot	
# List of bootstrapping results; only available when argument boot.val=TRUE

s_satpls$outer_model
s_satpls$inner_model
s_satpls$path_coefs
s_satpls$scores
s_satpls$effects
s_satpls$crossloadings
dim(s_satpls$crossloadings)

# time_data <- system.time(
#   s_satpls <- splspm(dfCompile, path_matrix, blocks, modes, scheme="path",
#                      scaled=F, penalization = "enet",
#                      cross_validate =F, lambda = 1, maxiter = 100)
# )
# saveRDS(object = s_satpls,file = "Rcode/RealDataAnalysis/LOS_msPLS_results_11182021.RDS")
# Select the results
# s_satpls <- readRDS("Rcode/RealDataAnalysis/LOS_msPLS_results_11182021.RDS")
mydata_run1 <- readRDS("Rcode/RealDataAnalysis/Final_LOS_results_run1.RDS")
s_satpls <- mydata_run1$model
feature_slection <- s_satpls$outer_model %>% 
  filter(weight!=0) 

length(grep(pattern = "GEID",x = feature_slection$name,value = T))
length(grep(pattern = "RRBS",x = feature_slection$name,value = T))
length(grep(pattern = "meta",x = feature_slection$name,value = T))

# Chuan's result
Qiu_tb <- openxlsx::read.xlsx(xlsxFile = "../../Data/MultiOmics136/Paper/BiomarkerTable_TableS2.xlsx")

# Get the annotation
# gene selection
GE_selection <- feature_slection %>%  filter(block== "EXPL_X1") %>% 
  left_join(x = .,y = GE_annot,by=c("name" = "GEID"))
# Methylation
Metyl_selection <- feature_slection %>%  filter(block== "EXPL_X2") %>% 
  left_join(x = .,y = RRBS_annot,by=c("name" = "RRBSID"))
# get for metabolomics 
metabo_select <- feature_slection %>% filter(block== "EXPL_X3") %>% 
  left_join(x = .,y = metabo_annot,by=c("name" = "metaID"))

Q.GE <- Qiu_tb %>% filter(Type == "GeneExpression") %>% select("Name")
table(GE_selection$Gene_symbol%in%Q.GE$Name)
# Gene expression feature selection
# GE_selection %>% filter(Gene_symbol%in%Q.GE$Name) %>% openxlsx::write.xlsx(x = .,file = "GE_FeatureSlect_LOS119Sample.xlsx")
# GE_selection %>% openxlsx::write.xlsx(x = .,file = "GE_FeatureSlect_LOS119Sample.xlsx")


############### msPLSDA ##############
# X1 = gene expression
# X2 = Methylome
# X3 = Metabolites

my_score_data <- data.frame(X1=s_satpls$scores[,1],
                            X2=s_satpls$scores[,2],
                            X3=s_satpls$scores[,3])

my_score_data <- cbind(my_score_data,dfCompile$BMD)
colnames(my_score_data)[4] <- "myGroup"
library(nnet) # multinomial regression
# my_score_data$myGroup <- relevel(my_score_data$myGroup, ref = "1")
training_model <- multinom(myGroup ~ X1 + X2 + X3, data = my_score_data)
training_outcom_truth <- as.factor(my_score_data$myGroup)
pp <- fitted(training_model)
# Get the label based on the fitted prob.
# max_prob_position <- apply(pp, 1, which.is.max)
max_prob_position <- ifelse(pp>0.5,1,0)
max_prob_position <- factor(max_prob_position)
performance.output.training.overall <- caret::confusionMatrix(max_prob_position,training_outcom_truth)


############### Plot components for omics data ##############
library("ggplot2")
# figure: gene expression vs methylalome
sc.p1 <- my_score_data %>% transform(myGroup = as.factor(myGroup)) %>% 
  mutate(myGroup = ifelse(myGroup==0,"Low BMD","High BMD")) %>% 
  rename(Group="myGroup") %>% 
ggplot(data =.,mapping = aes(x = X1,y=X2,color=Group))+
  geom_point()+labs(title = "Discriminant analysis on LOS multi-omics data",
                    subtitle = "First components of trasncriptome and methylome in 119 Samples",
                    x="Gene expression component 1",y="Methylation component 1")+
  theme_light()
# ggsave(plot = sc.p1,filename = paste0("Rcode/RealDataAnalysis/LOS_DA_Plot_Gene_Methy_",Sys.Date(),".pdf"),width = 6,height = 6)

# figure: gene expression vs metabolite
sc.p2 <- my_score_data %>% transform(myGroup = as.factor(myGroup)) %>% 
  mutate(myGroup = ifelse(myGroup==0,"Low BMD","High BMD")) %>% 
  rename(Group="myGroup") %>% 
  ggplot(data =.,mapping = aes(x = X1,y=X3,color=Group))+
  geom_point()+labs(title = "Discriminant analysis on LOS multi-omics data",
                    subtitle = "First components of trasncriptome and metabolome in 119 Samples",
                    x="Gene expression component 1",y="Metabolite component 1")+
  theme_light()
# ggsave(plot = sc.p2,filename = paste0("Rcode/RealDataAnalysis/LOS_DA_Plot_Gene_Metabo_",Sys.Date(),".pdf"),width = 6,height = 6)

# figure: metabolites vs methylation
sc.p3 <- my_score_data %>% transform(myGroup = as.factor(myGroup)) %>% 
  mutate(myGroup = ifelse(myGroup==0,"Low BMD","High BMD")) %>% 
  rename(Group="myGroup") %>% 
  ggplot(data =.,mapping = aes(x = X2,y=X3,color=Group))+
  geom_point()+labs(title = "Discriminant analysis on LOS multi-omics data",
                    subtitle = "First components of methylome and metabolome in 119 Samples",
                    x="Methylation component 1",y="Metabolite component 1")+
  theme_light()
# ggsave(plot = sc.p3,filename = paste0("Rcode/RealDataAnalysis/LOS_DA_Plot_Metyl_Metabo_",Sys.Date(),".pdf"),width = 6,height = 6)

# 3D plots for three omics
# install.packages("scatterplot3d") # Install
library("scatterplot3d") # load
my_score_data.f <- my_score_data %>% transform(myGroup = as.factor(myGroup)) %>% 
  mutate(myGroup = ifelse(myGroup==0,"Low BMD","High BMD")) %>% 
  rename("Gene Expression" = "X1", "Methylome" = "X2" , "Metabolite" = "X3") %>%
  mutate(Group = ifelse(myGroup=="Low BMD",1,2)) %>% 
  mutate(myGroup = as.factor(myGroup))

colors <- c("#999999", "#E69F00")
colors <- colors[as.numeric(my_score_data.f$Group)]
s3d <- scatterplot3d(my_score_data.f[,1:3], pch = 16, color = colors,
              grid=TRUE, box=T)
legend(legend = levels(my_score_data.f$myGroup),x = -2,y=2.5,
       col =  c("#999999", "#E69F00"), pch = 16,horiz = F)
# # Export the figure 
# ggsave(filename = "Rcode/RealDataAnalysis/3D_scatter_BMD_07222022.pdf")

# 3D plots
# library("plotly")
# my_score_data %>% transform(myGroup = as.factor(myGroup)) %>% 
#   mutate(myGroup = ifelse(myGroup==0,"Low BMD","High BMD")) %>% 
#   rename(Group="myGroup") %>% 
#   plot_ly(data = .,x=~X1,y=~X2,z=~X3) %>% 
#   add_markers(color = ~myGroup)
# 
# 
# p <- my_score_data %>% transform(myGroup = as.factor(myGroup)) %>% 
#   mutate(myGroup = ifelse(myGroup==0,"Low BMD","High BMD")) %>% 
#   rename(Group="myGroup") %>% 
#   my_score_data %>% 
#   plot_ly(., x=~X1, y =~X2 , z =~X3, color = ~myGroup, 
#              type="scatter3d", mode="markers")
# 
# p

########## Count markers in each components ##########
rm(list = ls())
load("Rcode/RealDataAnalysis/LOS_U19_pilot.RData")
rm(list=c("GeneExpression","RRBS","metabo","tmp.data"))
library("tidyverse")
# Chuan's result
Qiu_tb <- openxlsx::read.xlsx(xlsxFile = "../../Data/MultiOmics136/Paper/BiomarkerTable_TableS2.xlsx")

# Three components
table(Qiu_tb$Type,Qiu_tb$Componenet)
# msPLSDA:  the selection results
s_satpls <- readRDS("Rcode/RealDataAnalysis/LOS_msPLS_results_11182021.RDS")

feature_slection <- s_satpls$outer_model %>% 
  filter(weight!=0) 

length(grep(pattern = "GEID",x = feature_slection$name,value = T))
length(grep(pattern = "RRBS",x = feature_slection$name,value = T))
length(grep(pattern = "meta",x = feature_slection$name,value = T))

# Chuan's result
Qiu_tb <- openxlsx::read.xlsx(xlsxFile = "../../Data/MultiOmics136/Paper/BiomarkerTable_TableS2.xlsx")

# Get the annotation for msPLSDA
# gene selection
GE_selection <- feature_slection %>%  filter(block== "EXPL_X1") %>% 
  left_join(x = .,y = GE_annot,by=c("name" = "GEID"))
# Methylation
Metyl_selection <- feature_slection %>%  filter(block== "EXPL_X2") %>% 
  left_join(x = .,y = RRBS_annot,by=c("name" = "RRBSID")) %>% transform('CpG ID' = gsub(pattern = "_",replacement = "_Chr",x = `CpG ID`))
# get for metabolomics 
metabo_annot <- metabo_annot %>% transform("CNAME" = word(string = metabo_annot$Name,start = 1,end = 1,sep = "Δ"))
metabo_select <- feature_slection %>% filter(block== "EXPL_X3") %>% 
  left_join(x = .,y = metabo_annot,by=c("name" = "metaID"))
cat(metabo_select$CNAME,"\n")
for (i in 1:length(metabo_select$CNAME)) {
  # print(i)
  cat(metabo_select$CNAME[i],"\n")

}

#
# count.annot <- strsplit(x = metabo_annot$Name,split = "@") %>% lapply(.,FUN = length) %>% unlist 
# tmp.count <- which(count.annot == 1)
# uni.metabo_annot <- word(string = metabo_annot$Name[tmp.count],start = 1,end = 1,sep = "Δppm=")
# metabo_annot[tmp.count,"CNAME"] <- uni.metabo_annot



## check the overlap between msPLSDA and 
# install.packages("ggvenn")
library("ggvenn")
library("ggVennDiagram")
# https://www.datanovia.com/en/blog/venn-diagram-with-r-or-rstudio-a-million-ways/
GE_Qiu <- Qiu_tb %>% filter(Type == "GeneExpression") %>% select(Name) %>%  pull
Metyl_Qiu <- Qiu_tb %>% filter(Type == "Metylation") %>% select(Name) %>%  pull

# check multiple annotations
Metabo_Qiu <- Qiu_tb %>% filter(Type == "Metabolite") %>% select(Name) %>%  pull

for(i in 1:length(Metabo_Qiu)){
  # print(i)
  tmp.p <- grep(pattern = Metabo_Qiu[i],x = metabo_select$Name,fixed = T,value = T)
  tmp.pos <- grep(pattern = Metabo_Qiu[i],x = metabo_select$Name,fixed = T)
  if(!identical(tmp.p, character(0))){
    print(tmp.pos)
    cat(Metabo_Qiu[i],":")
    print(tmp.p)  
    cat("\n")
    metabo_select[tmp.pos,"Name"] <- Metabo_Qiu[i]
  }
}

# input for ggvenn is a list
x <- list(
  msplsda = c(GE_selection$Gene_symbol,Metyl_selection$CpG.ID,metabo_select$Name),
  smdcca = Qiu_tb$Name)
# ggvenn(data = x,stroke_size = 0.5,set_name_size = 4)
ggVennDiagram(x, label_alpha = 5)+
  ggplot2::scale_fill_gradient(low="lightblue",high = "yellow")

x_gene <- list(
  msplsda_Gene=GE_selection$Gene_symbol,
  smdcca_Gene=GE_Qiu
)
ggVennDiagram(x_gene, label_alpha = 5)+
  ggplot2::scale_fill_gradient(low="lightblue",high = "yellow")

x_metyl <- list(
  msplsda_CpG=Metyl_selection$CpG.ID,
  smdcca_CpG=Metyl_Qiu
)

ggVennDiagram(x_metyl, label_alpha = 5)+
  ggplot2::scale_fill_gradient(low="lightblue",high = "yellow")


x_metabo <- list(
  msplsda_metabolite=metabo_select$Name,
  smdcca_metabolite =Metabo_Qiu
)

ggVennDiagram(x_metabo, label_alpha = 5)+
  ggplot2::scale_fill_gradient(low="lightblue",high = "yellow")

x2 <- list(
  msplsda_Gene=GE_selection$Gene_symbol,
  msplsda_CpG=Metyl_selection$CpG.ID,
  msplsda_metabolite=metabo_select$Name,
  smdcca_Gene=GE_Qiu,
  smdcca_CpG=Metyl_Qiu,
  smdcca_metabolite =Metabo_Qiu
)
ggVennDiagram(x = x2)+
  ggplot2::scale_fill_gradient(low="lightblue",high = "yellow")

### Pathway analysis ###
GWAS_bone <-openxlsx::read.xlsx(xlsxFile = "/Users/raydai/Google\ Drive/Research/MultiOmicsFA/EnetPLS_project/EnetPLS/Rcode/RealDataAnalysis/Enrichr/gwas_catalog_v1.0.2-asso_e105_r2021-12-21_bone_selected.xlsx")
GE_selection$Gene_symbol
GWAS_bone_gene <- GWAS_bone$REPORTED.GENE.S.[which(!is.na(GWAS_bone$REPORTED.GENE.S.))]
x_gwas <- list(GE_msplasda = GE_selection$Gene_symbol,
               smdcca_Gene=GE_Qiu,
               GWAS_gene = GWAS_bone_gene)

ggVennDiagram(x = x_gwas)+
  ggplot2::scale_fill_gradient(low="lightblue",high = "yellow")




###################### component 2 ######################
mydata_run1 <- readRDS("Rcode/RealDataAnalysis/Final_LOS_results_run2.RDS")
s_satpls <- mydata_run1$model

my_score_data <- data.frame(X1=s_satpls$scores[,1],
                            X2=s_satpls$scores[,2],
                            X3=s_satpls$scores[,3])

my_score_data <- cbind(my_score_data,dfCompile$BMD)
colnames(my_score_data)[4] <- "myGroup"
library("ggplot2")
# figure: gene expression vs methylalome
sc.p1 <- my_score_data %>% transform(myGroup = as.factor(myGroup)) %>% 
  mutate(myGroup = ifelse(myGroup==0,"Low BMD","High BMD")) %>% 
  rename(Group="myGroup") %>% 
  ggplot(data =.,mapping = aes(x = X1,y=X2,color=Group))+
  geom_point()+labs(title = "Discriminant analysis on LOS multi-omics data",
                    subtitle = "Second components of trasncriptome and methylome in 119 Samples",
                    x="Gene expression component 2",y="Methylation component 2")+
  theme_light()
sc.p1
# ggsave(plot = sc.p1,filename = paste0("Rcode/RealDataAnalysis/LOS_DA_Plot_Gene_Methy_",Sys.Date(),"_comp2.pdf"),width = 6,height = 6)

# figure: gene expression vs metabolite
sc.p2 <- my_score_data %>% transform(myGroup = as.factor(myGroup)) %>% 
  mutate(myGroup = ifelse(myGroup==0,"Low BMD","High BMD")) %>% 
  rename(Group="myGroup") %>% 
  ggplot(data =.,mapping = aes(x = X1,y=X3,color=Group))+
  geom_point()+labs(title = "Discriminant analysis on LOS multi-omics data",
                    subtitle = "Second components of trasncriptome and metabolome in 119 Samples",
                    x="Gene expression component 2",y="Metabolite component 2")+
  theme_light()
# ggsave(plot = sc.p2,filename = paste0("Rcode/RealDataAnalysis/LOS_DA_Plot_Gene_Metabo_",Sys.Date(),"_comp2.pdf"),width = 6,height = 6)

# figure: metabolites vs methylation
sc.p3 <- my_score_data %>% transform(myGroup = as.factor(myGroup)) %>% 
  mutate(myGroup = ifelse(myGroup==0,"Low BMD","High BMD")) %>% 
  rename(Group="myGroup") %>% 
  ggplot(data =.,mapping = aes(x = X2,y=X3,color=Group))+
  geom_point()+labs(title = "Discriminant analysis on LOS multi-omics data",
                    subtitle = "Second components of methylome and metabolome in 119 Samples",
                    x="Methylation component 2",y="Metabolite component 2")+
  theme_light()
# ggsave(plot = sc.p3,filename = paste0("Rcode/RealDataAnalysis/LOS_DA_Plot_Metyl_Metabo_",Sys.Date(),"_comp2.pdf"),width = 6,height = 6)

