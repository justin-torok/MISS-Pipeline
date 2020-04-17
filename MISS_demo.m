% This script walks through the data processing steps necessary to generate
% the processed scRNAseq and ISH data stored in PresetInputs.mat, which are
% then used to perform cell count inference and parameter estimation. The
% resulting information is stored in a variable called "outstruct" with the
% index of maximum performance indicated by "idx". These are necessary
% inputs to many of the figure generation functions.

%% Grouping raw exon and intron counts into a hierarchical data struct per AIBS scRNAseq cell type labels
classstruct = scRNAseq_Data_Extract();

%% Generating the mean scRNAseq expression per cell type for consensus genes
[genevct, classkey, gene_names] = Cell_Type_Data_Extract(classstruct);
% Stored in PresetInputs.mat

%% Generating the per-voxel ISH gene expression for consensus genes
regvgene = ISH_Data_Extract(classstruct);
% Stored in PresetInputs.mat

%% Performing cell type inference and parameter optimization
directory = [cd filesep 'MatFiles'];
load([directory filesep 'PresetInputs.mat']);
clearvars -except regvgene genevct classkey gene_names
method = 'MRx3';
testnG = 300:120:900;
if strcmp('MRx3',method)
    testlambda = [150,500];
    [outstruct,idx] = lambda_ParameterFitter(regvgene, genevct, gene_names, testnG, testlambda);
else
    [outstruct,idx] = nG_ParameterFitter(regvgene, genevct, gene_names, method, testnG);
end
% Saved in MISS_demo_output.mat with preset parameter values above

%% A sampling of results from the MISS manuscript
directory = [cd filesep 'MatFiles' filesep 'MISS'];
load([directory filesep 'MISS_demo_output.mat']); 
savenclose = 0;
Figure_2c_bigscatter(outstruct,idx,savenclose,directory);
Figure_3b_interneuron(outstruct,idx,savenclose,directory);
Figure_3de_glia(outstruct,idx,savenclose,directory);
Figure_4ab_taulayerslice(outstruct,idx,[25,31,36],savenclose,directory);
Figure_5ab_exinhplots(outstruct,idx,savenclose);
randstruct = Rand_Index_Calc(outstruct,idx,'fore',directory);
Figure_6bc_aristats(randstruct,savenclose);

%% Plotting voxel renderings of cell type distributions
typeinds = [3,6,7,10,15]; % Refer to classkey
Cell_Type_Brainframe(outstruct,idx,typeinds,savenclose,[],directory);

