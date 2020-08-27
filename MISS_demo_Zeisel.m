% This script walks through the data processing steps necessary to generate
% the processed scRNAseq and ISH data stored in PresetInputs.mat, which are
% then used to perform cell count inference and parameter estimation. The
% resulting information is stored in a variable called "outstruct" with the
% index of maximum performance indicated by "idx". These are necessary
% inputs to many of the figure generation functions.

%% Grouping raw exon and intron counts into a hierarchical data struct per AIBS scRNAseq cell type labels
matdir = '/Users/justintorok/Documents/MATLAB/CellType/MatFiles/MISS'; % indicate folder where mat files are located
% classstruct = scRNAseq_Data_Extract(matdir);

%% Generating the mean scRNAseq expression per cell type for consensus genes
% [genevct] = Cell_Type_Data_Extract(classstruct,matdir);
% [C_indivcells, classkey, gene_names, ct_labvec] = Cell_Type_Extract_IndivCells(classstruct,matdir);
% Stored in PresetInputs.mat

%% Generating the per-voxel ISH gene expression for consensus genes
% regvgene = ISH_Data_Extract(classstruct,matdir);
regvgene = ISH_Data_Extract_Zeisel(matdir);
% Stored in PresetInputs.mat

%% Performing cell type inference and parameter optimization
load([matdir filesep 'Zeisel_extract.mat'],'C_indivcells','classkey','ct_labvec','entrez_names','meanexprmat');
% clearvars -except regvgene genevct classkey gene_names ct_labvec C_indivcells matdir
method = 'MRx3';
k = 5;
testlambda = 250;
testnG = 1168;
genevct = meanexprmat.'; C_indivcells = C_indivcells.'; ct_labvec = ct_labvec.';
sig = 4400;
% [~,~,~,~,mrx3inds] = GeneSelector(genevct,regvgene,entrez_names,3803,testlambda,method);
% parpool(4);
[fitstruct, outstruct] = nG_ParameterFitter_Zeisel(regvgene, genevct, entrez_names, method, testnG, testlambda, k, C_indivcells, ct_labvec, sig, matdir);
% delete(gcp('nocreate'))
% Saved in MISS_demo_output.mat with preset parameter values above
save([matdir filesep 'tasic_colAMD.mat'],'outstruct','fitstruct');

%% A sampling of results from the MISS manuscript
% directory = [cd filesep 'MatFiles' filesep 'MISS'];
% load([directory filesep 'MISS_demo_output.mat']); 
% savenclose = 0;
% Figure_2c_bigscatter(outstruct,idx,savenclose,directory);
% Figure_3b_interneuron(outstruct,idx,savenclose,directory);
% Figure_3de_glia(outstruct,idx,savenclose,directory);
% Figure_4ab_taulayerslice(outstruct,idx,[25,31,36],savenclose,directory);
% Figure_5ab_exinhplots(outstruct,idx,savenclose);
% randstruct = Rand_Index_Calc(outstruct,idx,'fore',directory);
% Figure_6bc_aristats(randstruct,savenclose);

%% Plotting voxel renderings of cell type distributions
% typeinds = [3,6,7,10,15]; % Refer to classkey
% Cell_Type_Brainframe(outstruct,idx,typeinds,savenclose,[],directory);

