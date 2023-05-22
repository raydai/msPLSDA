# Metadata ####
# Working example and simulation study for msPLS
#
#
# The R implementation of msPLS relies on the PLS-PM: Partial Least Squares Path Modeling R package from Sanchez et al. (2009).
# The package is modified so that it fits the requirements of big data analyis, as described in our manuscript.
# First, we provide an example of how msPLS can be applied to 3 datasets to find the relevant variables, and we also provive a section
# describing one of the data analysis studies we used in the manuscript. 

# Load the functions
library(R.utils) # Utility functions useful when programming and developing R packages.
library(mvtnorm)
library(tester) # tester allows you to test characteristics of common R objects.
library(turner) # Package designed for working with vectors and lists of vectors, mainly for turning them into other indexed data structures.
library(diagram) # Functions for Visualising Simple Graphs (Networks), Plotting Flow Diagrams
library(shape) # Functions for plotting graphical shapes such as ellipses, circles, cylinders, arrows, ...
library(amap)  # Tools for Clustering and Principal Component Analysis (With robust methods, and parallelized functions).
sourceDirectory("2018_msPLS-master/Example/sPLSPM/R/")
files.sources = list.files(path = "2018_msPLS-master/Example/sPLSPM/R/",full.names = T)
sapply(files.sources, source)
# an introductory example ####

#function for multiset multivariable data simulation
simulate_data <- function(sample_size = 100, 
                          nr_variables = c(50,50,20), 
                          nr_correlated_vars = c(5,5,5),
                          weight_correlations = c(0.7,0.5,0.3),
                          thetas = c(0.7,0.8,0.5),
                          sigma_correlation = 0.3){
  
  N = sample_size
  
  p1=nr_variables[1]
  p2=nr_variables[2]  
  p3=nr_variables[3]
  k=nr_correlated_vars[1]   # number of correlated varaibles
  
  w1_cor = weight_correlations[1]
  w2_cor = weight_correlations[2]
  w3_cor = weight_correlations[3]
  
  theta1 = thetas[1]
  theta2 = thetas[2]
  
  Sigma = diag((p1+p2))
  Sigma[Sigma == 0] = 0   # ?? What is this??
  #cor(X) 
  # Sigma = diag(10)
  # Sigma[(1:5)+5,(1:5)] = Sigma[(1:5),(1:5)+5] = sigma_correlation
  #First k correlated with the first value from Y.
  Sigma[(1:k)+p1,(1:k)] = Sigma[(1:k),(1:k)+p1] = sigma_correlation
  
  data = rmvnorm(N, ,Sigma)
  
  X1 = data[,1:p1]
  X2 = data[,(p1+1):(p1+p2)]
  
  # corrplot::corrplot(cor(X1[,1:10], X2[,1:10]))
  # corrplot::corrplot(cor(X1, X2))
  
  # Related features in X1
  w1 = rep(0,p1)
  w1[1:k] = rep(w1_cor,k)
  # Related features in X2
  w2 = rep(0,p2)
  w2[1:k] = rep(w2_cor,k)
  
  ksi1 = X1%*%w1     # [100, 50] X [50 X 1]
  ksi2 = X2%*%w2
  
  # solve(t(X1)%*%X1) %*% t(X1) %*% ksi1
  
  meanx = (theta1*ksi1) + (theta2*ksi2)
  sdx = (theta1)^2 + (theta2)^2
  
  ksi3 = rnorm(N,meanx,sqrt(abs(1-sdx)))
  
  w3 = rep(0,p3)
  w3[1:k] = rep(w3_cor,k)
  
  X3 = matrix(rep(0,N*p3),N,p3)
  
  for (j in 1:k) {
    sdx=sqrt(1-(w3)^2)
    X3[,j]=rnorm(N,meanx,sdx)      # sample values for the manifest variables
  }
  
  for (j in (k+1):p3) {
    X3[,j]=rnorm(N,0,1)      # sample values for the manifest variables
  }
  
  data$X1 <- X1
  data$X2 <- X2
  data$X3 <- X3
  
  return( data )
}

sample_size = 100 
nr_variables = c(50,50,50)
nr_correlated_vars = c(10)
weight_correlations = c(0.7,0.6,0.3)
thetas = c(0.8,0.7,0.3)
sigma_correlation = 0.5

# for replication
nr_seed = 73988636
set.seed(nr_seed)


#generate 3 datasets each with 10 relevant variables placed at the first variables' positions
data <- simulate_data(sample_size = sample_size, 
                      nr_variables = nr_variables, 
                      nr_correlated_vars = nr_correlated_vars,
                      weight_correlations = weight_correlations,
                      thetas = thetas,
                      sigma_correlation = sigma_correlation)
X <- data$X1
Y <- data$X2
Z <- data$X3

Data <- cbind(X,Y,Z)

dim(X)
dim(Y)
dim(Z)

dim(Data)

EXPL_X = c(0,1,0)
EXPL_Y = c(1,0,0)
RESP_Z = c(1,1,0)
path_matrix = rbind(EXPL_X, EXPL_Y,RESP_Z)


# blocks of outer model
blocks = list(1:dim(X)[2], 
              dim(X)[2]+1:dim(Y)[2],
              dim(X)[2]+dim(Y)[2]+1:dim(Z)[2])

dim(Data[,blocks[[1]]])
modes = c("B","B","A")


time_data <- system.time(
    s_satpls <- splspm(Data, path_matrix, blocks, modes, scheme="path",
                       scaled=T, penalization = "enet", nonzero = c(5,10,15),
                       lambda = c(0.1,1), maxiter = 100, cross_validate = T, nr_subsets = 10)
)

summary(s_satpls)

#the first 10 variables are assigned with non-zero weights from dataset X
s_satpls$outer_model[s_satpls$outer_model[,2]=="EXPL_X",]$weight

# and the first 10 varaibes are assigned with non-zero weights from dataset Y
s_satpls$outer_model[s_satpls$outer_model[,2]=="EXPL_Y",]$weight

# the regression weights for dataset Z are not penalized, but the first 10 variables are 
# assigned with higher regression weights then the rest
s_satpls$outer_model[s_satpls$outer_model[,2]=="RESP_Z",]$weight


################################################
# splspm example (Ray)
################################################
Data=Data
path_matrix
blocks
modes
scheme="path"
scaled=T
penalization = "enet"
nonzero = c(5,10,15)
lambda = c(0.1,1)
maxiter = 100
cross_validate = T
nr_subsets = 10

tol = 0.000001
scaling = NULL

plscomp = NULL
#' @param plscomp optional vector indicating the number of PLS components
#' (for each block) to be used when handling non-metric data
#' (only used if \code{scaling} is provided)
#' 
#' @param br number bootstrap resamples. Used only
#' when \code{boot.val=TRUE}. When \code{boot.val=TRUE}, the default number of
#' re-samples is 100.
boot.val = FALSE
br = NULL
#' @return \item{data}{Data matrix containing the manifest variables used in the
#' model. Only available when \code{dataset=TRUE}}
dataset = TRUE
alpha = 0.5
warning_non_convergence = TRUE





