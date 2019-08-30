# Description
These programs are made for "Tropical Statistics for Phylogenetic Trees" paper.

1. func_ssh.R: The defined functions used in plot.R and r_square.R.
2. plot.R:     A script used to generate Figure 2 and Figure 3 in the paper.
3. r_square.R: A function used to calculate r square in tropical space.

Two principal components analysis (PCA) in tree space mentioned in the paper:

* Tropical PCA: https://arxiv.org/abs/1710.02682.
* BHV PCA: https://arxiv.org/abs/1609.03045.

# Required R package

* ape
* phangorn
* parallel
* lpSolveAPI

# Environment
We use parallel computing because the method used in the paper requires a lot of calculations. Only a Windows version is available right now. Users can modify the source code to create the Linux version according to their needs.
