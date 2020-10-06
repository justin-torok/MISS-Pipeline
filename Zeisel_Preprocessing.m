% This script walks through the data preprocessing steps necessary to
% create the input data matrices for the MISS pipeline using the Zeisel et
% al. dataset (Zeisel et al., 2018 - http://mousebrain.org/downloads.html).
% This requires that, in addition to the other dependencies in the provided
% 'MatFiles' folder, the user download the 'L5_All.agg.loom' and
% 'L5_All.loom' files from the website above. It also uses the AIBS mouse
% gene expression atlas (Lein et al., 2007) and removes all nonoverlapping
% genes between the two datasets.

%% Generating the mean scRNAseq expression per cell type for consensus genes
matdir = '/Users/justintorok/Documents/MATLAB/CellType/MatFiles/MISS'; % indicate folder where mat files are located
[genevct, C_indivcells, classkey, gene_names, ct_labvec] = Cell_Type_Data_Extract_Zeisel(matdir);

%% Generating the per-voxel ISH gene expression for consensus genes
voxvgene = ISH_Data_Extract_Zeisel(matdir);

% All outputs stored in Zeisel_Inputs.mat, along with other useful 
% datafiles related to the regional parcellation of the ISH atlas.