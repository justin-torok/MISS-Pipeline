function [fitstruct, outstruct] = nG_ParameterFitter_Zeisel(voxvgene_, genevct_, gene_names_, method_, ng_param_list_, lambda_, k_, C_indivcells_, ct_labvec_, sigmas_, matdir_)
% This function performs a brute-force parameter sweep to find the maximal
% performance of a given subset selection method. The user MUST supply
% voxvgene and genevct, which are outputs of ProcessedData_Generator.m, and
% method, a string identifying a specific subset selection method. The
% outputs of this function are outstruct, which contains all of the
% inference and fitting-related information necessary to recreate plots in
% Mezias et al, 2020, and peakind, which identifies which index of
% outstruct produced maximal sumfit. 
%
% This function calls upon various dependent functions that perform
% different aspects of the data processing and quantitative assessment. For
% MRx3, lambda is fit in a separate function (lambda_ParameterFitter.m)
% that calls upon this function.

% Default nG parameter lists
if nargin < 11
    matdir_ = cd;
    if nargin < 10
        sigmas_ = 4400;
        if nargin < 9
            ct_labvec_ = [];
            if nargin < 8
                C_indivcells_ = [];
                if nargin < 7
                    k_ = 4;
                    if nargin < 6
                        lambda_ = 150;
                        if nargin < 5
                            if nargin < 4
                                error(sprintf('Error. \nUser must supply a voxel x gene matrix, a gene x cell type matrix, a gene names cell array, and a subsetting method'))
                            end
                            if ismember(method_, {'MRx3','mRMR','colAMD'})
                                ng_param_list_ = [100:70:1500, 1750:250:3500, 3855]; 
                            elseif strcmp(method_, 'DBSCAN')
                                ng_param_list_ = [0.0001, 0.0004:0.0001:0.001, 0.002:0.001:0.005, 0.015:0.01:0.195]; 
                            elseif strcmp(method_, 'Entropy')
                                ng_param_list_ = [0.25,0.5:0.125:2.75,2.8:0.05:3.25];
                            elseif strcmp(method_, 'Zeisel')
                                ng_param_list_ = [2:5, 6:2:50]; % made up something here
                            else
                                error('Error. \n%s is an incorrect subsetting method identifier',method_)
                            end
                        end    
                    end
                end
            end
        end
    end
end

fitstruct = struct;
outstruct = struct;

fprintf('Initializing preloaded gene indices\n')
if strcmp(method_,'MRx3')
    preloadinds = MRx3_Selector(genevct_,voxvgene_,max(ng_param_list_),lambda_);
else
    preloadinds = [];
end

for i = 1:length(ng_param_list_)
    fprintf('Determining subset, nG parameter value %d/%d\n',i,length(ng_param_list_))
    param = ng_param_list_(i);
    % Create reduced versions of voxvgene and genevct according to the
    % method and parameter specified by the user
    [~,~,nGen,~,C_ind_red] = GeneSelector_IndivCells(genevct_,voxvgene_,C_indivcells_,gene_names_,param,lambda_,method_,preloadinds);
    
    % Calculate classification error
    tic
    fprintf('Determining GMM classification error, nG parameter value %d/%d\n',i,length(ng_param_list_))
    savegroups = 1;
    crossval_ = 0;
    gmmstruct = GMM_Nearest_Neighbor_Posterior(C_ind_red, ct_labvec_, k_, crossval_, savegroups);
    fitstruct(i).gmmstruct = gmmstruct;
    fitstruct(i).lambda = lambda_;
    fitstruct(i).crossval = crossval_;
    fitstruct(i).nGen = nGen;
    fitstruct(i).sigmas = sigmas_;
    toc
    fprintf('Done, GMM fitting, nG parameter value %d/%d\n',i,length(ng_param_list_))
end

% Elbow determination for nG range supplied
fprintf('Determining optimal nG value\n');
neglogpriors = zeros(length(fitstruct),length(sigmas_));
negloglikelihoods = neglogpriors;
for i = 1:length(fitstruct)
    for j = 1:length(sigmas_)
        negloglikelihoods(i,j) = -(fitstruct(i).gmmstruct.gmmpost^2);
        neglogpriors(i,j) = (ng_param_list_(i)^2)/(2*sigmas_(j)^2);
    end
    fitstruct(i).negloglikelihood = negloglikelihoods(i,1);
    fitstruct(i).neglogpriors = neglogpriors(i,:);
end

neglogposteriors = negloglikelihoods + neglogpriors;
[~,mininds] = min(neglogposteriors);
nG_opts = ng_param_list_(mininds);

for i = 1:length(nG_opts)
    outstruct(i).nGen = nG_opts(i);
    % Infer cell density per voxel in arbitrary units
    tic
    fprintf('Nonnegative matrix inversion, optimal nG parameter value %d/%d\n',i,length(nG_opts))
    [E_red,C_red] = GeneSelector(genevct_,voxvgene_,gene_names_,nG_opts(i),lambda_,method_,preloadinds);
    B = CellDensityInference(E_red,C_red);
    toc
    %     outstruct(i).Bvals = B; 
    if strcmp(method_,'MRx3')
        outstruct(i).lambda = lambda_;
    end

    % Convert arbitrary densities to counts per voxel
    Bcorrected = Density_to_Counts(B,matdir_);
    outstruct(i).corrB = Bcorrected;

    % Sum and average over CCF regions
    [sumB,meanB] = Voxel_To_Region(Bcorrected,matdir_);
    outstruct(i).Bsums = sumB; % total cells per region
    outstruct(i).Bmeans = meanB; % mean cell count per region
end
end
