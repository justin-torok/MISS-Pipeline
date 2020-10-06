% This script walks through the data processing steps necessary to generate
% the processed scRNAseq and ISH data stored in PresetInputs.mat, which are
% then used to perform cell count inference and parameter estimation. The
% resulting information is stored in a variable called "outstruct" with the
% index of maximum performance indicated by "idx". These are necessary
% inputs to many of the figure generation functions.

%% Performing hyperparameter optimization
matdir = '/Users/justintorok/Documents/MATLAB/CellType/MatFiles/MISS'; % indicate folder where mat files are located
% method = 'MRx3';
% k = 5;
% cv = 1;
% nG_vals = [50:50:300,301:1200,1255:50:2155,2255:100:3655];
% lambda_list = 0:50:350;
% sigma_list = 3000:100:6000;
% [lambda_opt,sigma_opt,ngstar_opt] = Hyperparameter_Optimizer(method,...
%                                                 nG_vals, lambda_list,...
%                                                 sigma_list, k, cv, matdir);
                                            
% Outputs: lambda_opt = 250; sigma_opt = 4400; ngstar_opt = 529.

%% Create optimal cell type maps for the Tasic et al. dataset
% load([matdir filesep 'Tasic_Inputs.mat'],'voxvgene','genevct','gene_names',...
%     'C_indivcells','ct_labvec');
% lambda_opt = 250; 
% sigma_opt = 4400; 
% ngstar_opt = 529;
% method = 'MRx3';
% k = 5;
% cv = 1;
% [fitstruct, outstruct] = nG_ParameterFitter(voxvgene, genevct,...
%                                     gene_names, method, C_indivcells,...
%                                     ct_labvec, ngstar_opt, lambda_opt,...
%                                     sigma_opt, k, cv, matdir);

% Stored in the file 'tasic_l250_ng529.mat'

%% Create optimal cell type maps for the Zeisel et al. dataset
% load([matdir filesep 'Tasic_Inputs.mat'],'voxvgene','genevct','gene_names',...
%     'C_indivcells','ct_labvec');
% lambda_opt = 250; 
% sigma_opt = 4400; 
% method = 'MRx3';
% k = 5;
% cv = 1;
% nG_vals = [50:50:300,301:1200,1255:50:2155,2255:100:3655];
% [fitstruct, outstruct] = nG_ParameterFitter(voxvgene, genevct,...
%                                     gene_names, method, C_indivcells,...
%                                     ct_labvec, nG_vals, lambda_opt,...
%                                     sigma_opt, k, cv, matdir);

% Stored in the file 'zeisel_l250_ng1168.mat'

%% A sampling of results from the MISS manuscript
load([matdir filesep 'tasic_l250_ng529.mat'],'outstruct');
savenclose = 0;
Figure_2ab_posteriors_hyperparam(savenclose,matdir);
Figure_3ab_methodcomp(savenclose,matdir);
Figure_4b_interneuron(outstruct,1,savenclose,matdir);
% Figure_3de_glia(outstruct,idx,savenclose,directory);
% Figure_4ab_taulayerslice(outstruct,idx,[25,31,36],savenclose,directory);
% Figure_5ab_exinhplots(outstruct,idx,savenclose);
% randstruct = Rand_Index_Calc(outstruct,idx,'fore',directory);
% Figure_6bc_aristats(randstruct,savenclose);

%% Plotting voxel renderings of cell type distributions
typeinds = [3,6,7,10,15]; % Refer to classkey
Cell_Type_Brainframe(outstruct,idx,typeinds,savenclose,[],directory);

