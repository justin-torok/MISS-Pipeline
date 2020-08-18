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
% Stored in PresetInputs.mat

%% Performing cell type inference and parameter optimization
load([matdir filesep 'PresetInputs.mat']);
clearvars -except regvgene genevct classkey gene_names ct_labvec C_indivcells matdir
method = 'MRx3';
% testnG = [50:50:450,460:10:1600,1650:50:2500,2600:100:3700,3855];
% testnG = [50:50:300,301:600,610:10:1200];
testnG = [120:100:820];
% testnG = [369,383,391,425];
% costfun = 'SumFit';
k = 5;
testlambda = 150;
sig = [3800,4400,5000];
[fitstruct, outstruct] = nG_ParameterFitter_Zeisel(regvgene, genevct, gene_names, method, testnG, testlambda, k, C_indivcells, ct_labvec, sig, matdir);


% if strcmp('MRx3',method)
%     testlambda = 150;
%     [outstruct,idx] = lambda_ParameterFitter(regvgene, genevct, gene_names, testnG, testlambda, costfun, k, C_indivcells, ct_labvec, matdir);
% else
%     testlambda = 150;
%     [outstruct,idx] = nG_ParameterFitter(regvgene, genevct, gene_names, method, testnG, testlambda, costfun, k, C_indivcells, ct_labvec, matdir);
% end
% save('/Users/justintorok/Documents/MATLAB/CellType/MatFiles/lambdatest_lda_150_10fold_1.mat', 'outstruct','idx','-v7.3');
% Saved in MISS_demo_output.mat with preset parameter values above

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

