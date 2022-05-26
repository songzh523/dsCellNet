# What is tsCellNet?
We present a cell-cell interaction networks inference method and its application for time-series RNA-seq data (termed tsCellNet). 
It identifies interactions according to protein localizations and reinforces them via dynamic time warping within each time point and between adjacent time points, respectively. Then, fuzzy clustering of those interactions helps us refine key time points and connections. Compared to other published methods, our methods display high accuracy and high tolerance based on both simulated data and real data. Moreover, our methods can help identify the most active cell type and genes, which may provide a powerful tool to uncover the mechanisms underlying the organismal development and disease. 

# How to use tsCellNet?
## Step1:Installing tsCellNet
### Depends:
* R (R>=4.0.3)
* igraph(>=1.2.6); homologene(>=1.4.68); Hmisc(>=4.5.0); ggpubr(>=0.4.0); pheatmap(>=1.0.12); reshape2(>=1.4.4); ggalluvial(>=0.12.3); tseries(>=0.10.48); magrittr(>=2.0.1);

## Step2:Running tsCellNet
### Parameters:
** datloc: path of the three/two input files;
** expnam: filenames of expression file(mandatory); 
** typnam: filenames of cell-type file(mandatory); 
** timnam: filenames of time-point file(optional);
** resloc: path of result files; 
** mainlab: prefix name of the result files;
** species: species name;
** controltime: the starting time-point name;
** treattime: if or not need to compare cell-cell interactions with starting time-point, if not set the values as "total"
### Example:



