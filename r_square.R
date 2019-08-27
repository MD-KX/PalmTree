# Script used to calculate R-Square
#
# Author: Qiwen Kang

# Settings ----------------------------------------------------------------
library(lpSolveAPI)
library(parallel)
library(ape)
library(phangorn)
library(parallel)
setwd("E:r//181025tropMCMC/")
source("./program/func_ssh.R")
# Define the number of clusters will be used.
cl <- makeCluster(8) 


# Function ----------------------------------------------------------------
#` trees_ori: the original trees data, it should have a multiphy format.
#` comb_set:  the index of combination. 
#` pcs:       number of principal components. 
#` outliers:  the index of outlier trees in trees_ori.

r_squre <- function(trees_ori, comb_set, pcs = 3,  outliers = NA){
  
  # Here we remove outliers from original tree data set (Outliers are detected using KDETREE method). 
  if(sum(is.na(outliers)) == 0){
    trees_ori <- trees_ori[-c(outliers)]
  }
  
  # Since some branches are negative length, we set them equal to 0.
  for(i in 1:length(trees_ori)){
    if(sum(trees_ori[[i]]$edge.length) > 0){
      trees_ori[[i]]$edge.length[trees_ori[[i]]$edge.length < 0] = 0
    }
  }
  
  n <- length(trees_ori[[1]]$tip.label) 
  to <- trees_ori[[1]]$tip.label 
  
  distVec_all <- distMat(trees_ori, tipOrder = to) 
  N <- length(distVec_all) 
  D_all <- matrix(unlist(distVec_all), ncol=N)
  
  new_base <- D_all[, comb_set]
  proj_points <- parLapply(cl, distVec_all, project_pi , D_s = new_base) 
  tropical_dist_vec <- mapply(tropical_dist, distVec_all, proj_points) 
  sum_dist <- sum(tropical_dist_vec) 
 
  r_proj_data <- matrix(unlist(proj_points), nrow = length(proj_points), byrow = T)
  
  apicom_fermet <- fermatweberdistance(r_proj_data)
  r_apicom <- apicom_fermet/(sum_dist + apicom_fermet)
  r_apicom
}

