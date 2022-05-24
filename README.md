# What is tsCellNet?
We present a cell-cell interaction networks inference method and its application for time-series RNA-seq data (termed tsCellNet). 
It identifies interactions according to protein localizations and reinforces them via dynamic time warping within each time point and between adjacent time points, respectively. Then, fuzzy clustering of those interactions helps us refine key time points and connections. Compared to other published methods, our methods display high accuracy and high tolerance based on both simulated data and real data. Moreover, our methods can help identify the most active cell type and genes, which may provide a powerful tool to uncover the mechanisms underlying the organismal development and disease. 

# How to use tsCellNet?
## Step1:Installing tsCellNet
*R (R>=4.0.3)*
igraph(>=1.2.6); homologene(>=1.4.68); Hmisc(>=4.5.0); ggpubr(>=1.2.6); pheatmap(>=1.2.6); reshape2(>=1.2.6); ggalluvial(>=1.2.6); tseries(>=1.2.6); magrittr(>=1.2.6);

## Step2:Running tsCellNet


