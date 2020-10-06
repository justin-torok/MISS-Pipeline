% This script walks through the data processing steps necessary to generate
% the processed scRNAseq and ISH data stored in PresetInputs.mat, which are
% then used to perform cell count inference and parameter estimation. The
% resulting information is stored in a variable called "outstruct" with the
% index of maximum performance indicated by "idx". These are necessary
% inputs to many of the figure generation functions.

%% Grouping raw exon and intron counts into a hierarchical data struct per AIBS scRNAseq cell type labels
matdir = '/wynton/home/rajlab/jtorok/MATLAB/MatFiles/MISS'; % indicate folder where mat files are located
outdir = '/wynton/home/rajlab/jtorok/MATLAB/MatFiles/OutputFiles';
% classstruct = scRNAseq_Data_Extract(matdir);

%% Generating the mean scRNAseq expression per cell type for consensus genes
% [genevct] = Cell_Type_Data_Extract(classstruct,matdir);
% [C_indivcells, classkey, gene_names, ct_labvec] = Cell_Type_Extract_IndivCells(classstruct,matdir);
% Stored in PresetInputs.mat

%% Generating the per-voxel ISH gene expression for consensus genes
% regvgene = ISH_Data_Extract(classstruct,matdir);
% Stored in PresetInputs.mat

%% Performing cell type inference and parameter optimization
load([matdir filesep 'PresetInputs.mat']);
clearvars -except regvgene genevct classkey gene_names ct_labvec C_indivcells matdir
testlambda = 250;
sig = 4400;
k = 5;
method = 'MRx3';
testnG = [50:50:300,301:1200,1255:50:2155,2255:100:3655];
[fitstruct, outstruct] = nG_ParameterFitter(regvgene, genevct, gene_names, method, testnG, testlambda, k, C_indivcells, ct_labvec, sig, matdir);
sigs = 3000:25:6000;
[sumfits,minngs] = SumFit_Optimizer(fitstruct,method,sigs,matdir);
save([matdir filesep sprintf('tasic_MRx3_l%d.mat',testlambda)],'outstruct','fitstruct','sumfits','minngs');

