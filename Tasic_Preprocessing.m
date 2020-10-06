% This script walks through the data preprocessing steps necessary to
% create the input data matrices for the MISS pipeline using the Tasic et
% al. dataset (combined scRNAseq datasets from Tasic et al., 2018 and AIBS,
% 2018 - https://portal.brain-map.org/atlases-and-data/rnaseq). It also
% uses the AIBS mouse gene expression atlas (Lein et al., 2007) and removes
% all nonoverlapping genes between the two datasets.

%% Grouping raw exon and intron counts into a hierarchical data struct per AIBS scRNAseq cell type labels
matdir = '/Users/justintorok/Documents/MATLAB/CellType/MatFiles/MISS'; % indicate folder where mat files are located
classstruct = scRNAseq_Data_Extract_Tasic(matdir);

%% Generating the mean scRNAseq expression per cell type for consensus genes
[genevct] = Cell_Type_Data_Extract_Tasic(classstruct,matdir); % C
[C_indivcells, classkey, gene_names, ct_labvec] = Cell_Type_Data_Extract_IndivCells_Tasic(classstruct,matdir);

%% Generating the per-voxel ISH gene expression for consensus genes
voxvgene = ISH_Data_Extract_Tasic(classstruct,matdir);

% All outputs with the exception of classstruct stored in Tasic_Inputs.mat,
% along with other useful datafiles related to the regional parcellation of
% the ISH atlas.
