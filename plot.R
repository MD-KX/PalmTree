# R script used for plotting

# Settings ----------------------------------------------------------------
library(ape)
library(phangorn)
library(parallel)
setwd("C:/Users/kangq/Desktop/program")
# "func_ssh.R" file includes all the functions we need for calculating with tropical PCA 
source("./func_ssh.R")
cl <- makeCluster(6) # Define the number of clusters


# Tropical trees ----------------------------------------------------------
# Author: Qiwen Kang

# For more details, please check https://github.com/QiwenKang/tropPCA.

## Tropical principal components trees.
trees_ori <- read.tree("./data/NYh3n2_HA_20000_4_1993.txt")

# comb_set is the set includes the index of tree combination detected by the tropical PCA method.
comb_set <- c(8645, 15739, 6802) 
# outliers is the set includes the index of trees detected by KDETREES method.
outliers <- c() 

# Here we remove outliers from original tree data set.
if(sum(is.na(outliers)) == 0){
  trees_ori <- trees_ori[-c(outliers)]
}
  
for(i in 1:length(trees_ori)){
  if(sum(trees_ori[[i]]$edge.length) > 0){
    trees_ori[[i]]$edge.length[trees_ori[[i]]$edge.length < 0] = 0
  }
}
  
# Collect some basic information of the trees.
n <- length(trees_ori[[1]]$tip.label) 
pcs <- 3 
to <- trees_ori[[1]]$tip.label 

distVec_all <- distMat(trees_ori, tipOrder = to) 
N <- length(distVec_all)  
D_all <- matrix(unlist(distVec_all), ncol = N) 


trees <- trees_ori[comb_set]

# mm <- read.csv("E:/r/181025tropMCMC/output/comb_list_NYh3n2_HA_20000_4_1993.txt.csv")[,-1]
# cc <- c()
# for(i in 1:5){
#   new_base <- D_all[,mm[,i]]
#   proj_points <- parLapply(cl, distVec_all, project_pi , D_s = new_base)
#   tropical_dist_vec <- mapply(tropical_dist, distVec_all, proj_points)
#   cc[i] <- sum(tropical_dist_vec)
# }
# mm[,3]

base <- distMat(trees, tipOrder = to)
D_base <- matrix(unlist(base), ncol = pcs)
adj_proj_points <- lapply(distVec_all, project_pi,D_s = D_base)
proj_trees_matrix <- lapply(adj_proj_points, make.matrix, n = length(to), tips = to)
proj_trees <- lapply(proj_trees_matrix, upgma)
  
type <- seq(length(proj_trees))
for(i in seq(length(proj_trees))){
  if(i==type[i]) {
    for(j in seq(i, length(proj_trees))) {
      if(RF.dist(proj_trees[[i]],proj_trees[[j]]) == 0) {
        type[j] = i
      }
    } 
  }
} # Classify tree topogolies 
  
table(type) # Check how many tree topogogies are identified
  
# trees[[1]], trees[[2]], trees[[3]] are three principal components
# Here we combine different tree topologies with three principal components.
trees[[4]] <- proj_trees[[1]]


# Keep the tree topology.
for(i in 1:6) trees[[i]]$edge.length <- NULL 

par.old <- par(mar = c(0.2,0.25,1.1,0.25), xpd = NA, font.main = 1) # Define the parameters of the figure
on.exit({
  par(par.old)
  layout(1)
}) 
layout(matrix(1:6, nc = 3, byrow = TRUE))

# Replace the leaves label with new labels
tip.orig <- c("Sample1", "Sample2", "Sample3","Sample4")
tip.abbrv <- c("1", "2", "3", "4")
for (i in seq_along(trees)){
  tord <- match(trees[[i]]$tip.label, tip.orig)
  trees[[i]]$tip.label <- tip.abbrv[tord]
} 
treelabs <- c("PC1", "PC2", "PC3", "Topology 1 (20000)")

# Plot the trees.              
for(i in seq_along(trees)) ape::plot.phylo(trees[[i]], main=treelabs[i], type="c", direction="downwards",srt=90,adj=0.5,label.offset=0.4,
                                             cex = 2,cex.main=2,  edge.width = 2, font = 2, font.main=2) # Plot the trees
 

## Tropical triangle.
new_base <- D_all[,comb_set]
DD_base <- t(new_base)
D_base <- normalize.ultrametrices(DD_base)
proj_points_plot <- lapply(adj_proj_points, polytope_iso, D = DD_base)
  
proj_plot_norm <- lapply(proj_points_plot, normalize.proj)
proj_2D_plot_m <- matrix(unlist(proj_plot_norm), nrow = N, ncol = 3, byrow = T)
  
type_ind <- as.matrix(table(type))
  
tree_num <- as.numeric(rownames(type_ind))
freq <- type
  
# Here we need to define the color we need
for(i in 1:length(tree_num)){
  freq[type == tree_num[i]] <- i
}
  
  
par(mfrow = c(1,1))
par(mar = c(2,2,0.2,0.2), xpd = NA, font.main = 1)
k <- ncol(D_base)
plot(D_base[1,],D_base[2,], xlab = "", ylab = "",cex.axis = 1.3,font = 2)
for(i in 1:(k - 1)){
  for(j in (i + 1):k){
    tseg1 <- tropical.geodesic.dim.2(D_base[,i], D_base[,j])
    tseg2 <- tropical.geodesic.dim.2(D_base[,i], D_base[,j], flag=1)
    if(tseg1[[2]] < tseg2[[2]]) tseg <- tseg1
    else tseg <- tseg2
    segments(tseg[[1]][1,1], tseg[[1]][2,1], tseg[[1]][1,2], tseg[[1]][2,2], col = 'black')
    segments(tseg[[1]][1,2], tseg[[1]][2,2], tseg[[1]][1,3], tseg[[1]][2,3], col = 'black')
  }
}
points(x = proj_2D_plot_m[,2], y = proj_2D_plot_m[,3], pch = 16, cex = 0.75, col = freq)
  

# BHV metric  ---------------------------------------------------------
# Author: Qiwen Kang, Grady Weyenberg

# For more details, please check https://github.com/grady/geophyttertools.
library(geophyttertools)
library(ape)
simplex2cartesian <- function(simplex){
  A <- cbind(x=c(1, 0.5) ,y=c(0, sqrt(3/4)))
  as.matrix(simplex) %*% Atropical.geodesic.dim.2
}


proj.lines <- readLines("E:/working/5/NYh3n2_HA_20000_5_2013.colc")
proj.lines <- proj.lines[2:(length(proj.lines) - 1)]
if(sum(is.na(outliers[[index]])) == 0){
  proj.lines <- proj.lines[-c(outliers[[index]])]
}

tree_list_ori <- strsplit(proj.lines, split= " ")
tree_list <- lapply(tree_list_ori, function(x) x[6])

trees_bhv <- list()
for(i in 1:length(tree_list)){
  trees_bhv[[i]] <- as.list(ape::read.tree(text=tree_list[[i]]))
}


type <- seq(length(trees_bhv))
for(i in seq(length(trees_bhv))){
  if(i==type[i]) {
    for(j in seq(i, length(trees_bhv))) {
      if(RF.dist(trees_bhv[[i]],trees_bhv[[j]])==0) {
        type[j]=i
      }
    } 
  }
}

table(type)

tree1 <- trees_bhv[[1]] 
tree2 <- trees_bhv[[3]] 
tree3 <- trees_bhv[[5]] 
tree4 <- trees_bhv[[18]] 
tree5 <- trees_bhv[[68]] 

bhv <- read.tree("pathOfBHVtrees")

fechet <- as.list(ape::read.tree(text="((Sample2:11.4776712,(Sample1:4.1489723,Sample4:10.2672235):0.4042340):0.6535183,(Sample5:11.4709264,Sample3:4.1867112):0.6535183);"))

tip.orig <- c("Sample1", "Sample2", "Sample3","Sample4","Sample5")
tip.abbrv <- c("1", "2", "3", "4", "5")
tord <- match(tree1$tip.label, tip.orig)
tree1$tip.label <- tip.abbrv[tord]
tree1$edge.length <- NULL

tord <- match(tree2$tip.label, tip.orig)
tree2$tip.label <- tip.abbrv[tord]
tree2$edge.length <- NULL

tord <- match(fechet$tip.label, tip.orig)
fechet$tip.label <- tip.abbrv[tord]
fechet$edge.length <- NULL

for(i in 1:3){
  tord <- match(bhv[[i]]$tip.label, tip.orig)
  bhv[[i]]$tip.label <- tip.abbrv[tord]
  bhv[[i]]$edge.length <- NULL
}


layout(matrix(1:9, nc = 3, byrow=TRUE))
par.old <- par(mar=c(0.2,0.25,1.1,0.25), xpd=NA, font.main=1)
on.exit({
  par(par.old)
  layout(1)
})
plot.phylo(bhv[[1]], main = "Tree 1", type="c", direction="d", cex = 2, cex.main=1.5, edge.width = 2, font = 2, 
           font.main=2, srt=90,adj=0.5,label.offset=0.3)
plot.phylo(bhv[[2]], main = "Tree 2", type="c", direction="d", cex = 2, cex.main=1.5, edge.width = 2, font = 2, 
           font.main=2, srt=90,adj=0.5,label.offset=0.3)
plot.phylo(bhv[[3]], main = "Tree 3", type="c", direction="d", cex = 2, cex.main=1.5, edge.width = 2, font = 2, 
           font.main=2, srt=90,adj=0.5,label.offset=0.3)

plot.phylo(tree1, main = "Topology 1 (10586)",type="c", direction="d", cex = 2, cex.main=1.5, edge.width = 2, font = 2, 
           font.main=2, srt=90,adj=0.5,label.offset=0.3)
plot.phylo(tree2, main = "Topology 2 (4794)",type="c", direction="d", cex = 2, cex.main=1.5, edge.width = 2, font = 2, 
           font.main=2, srt=90,adj=0.5,label.offset=0.3)
plot.phylo(tree3, main = "Topology 3 (4616)",type="c", direction="d", cex = 2, cex.main=1.5, edge.width = 2, font = 2, 
           font.main=2, srt=90,adj=0.5,label.offset=0.3)
plot.phylo(tree4, main = "Topology 4 (3755)",type="c", direction="d", cex = 2, cex.main=1.5, edge.width = 2, font = 2, 
           font.main=2, srt=90,adj=0.5,label.offset=0.3)
plot.phylo(tree5, main = "Topology 5 (363)",type="c", direction="d", cex = 2, cex.main=1.5, edge.width = 2, font = 2, 
           font.main=2, srt=90,adj=0.5,label.offset=0.3)

plot.phylo(fechet, main = "Fréchet mean (10586)",type="c", direction="d", cex = 2, cex.main=1.5, edge.width = 2, font = 2, 
           font.main=2, srt=90,adj=0.5,label.offset=0.3)

## BHV triangle 
library(ggplot2)
# The ".trop" file is define by user. Basically, it is the output of the JAVA program, geophytterplus.DecomposeLFMTriangle.
background <- read.topologies("pathOfBHVoutputs")
# The ".colc" file is define by user. It is a part of the output of the JAVA program, geophytterplus.FitLFMTriangle.
points <- read.projections("pathOfBHVoutputs")
pp <- ggplot(background) + corner_labels(offset=0.03, size = 6) +
  guides(fill=guide_legend(title="Topology")) +
  #scale_fill_brewer(palette="Set3") +
  scale_fill_manual(values=c("#8dd3c7","#ffffb3","#bebada","#fb8072","#b3de69","#fdb462","#80b1d3","#fccde5","#d9d9d9")) +
  geom_point(aes(x=x,y=y), points, inherit.aes = FALSE, size=0.5) +
  theme(legend.position = c(.82,.76),legend.text=element_text(size=16, face = "bold"), legend.title = element_text(size = 16,face = "bold"))
pp
