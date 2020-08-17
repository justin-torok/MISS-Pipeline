function [outstruct, optind] = nG_ParameterFitter(voxvgene_, genevct_,...
                                    gene_names_, method_, ng_param_list_,...
                                    lambda_, costfun_, k_, C_indivcells_,...
                                    ct_labvec_, matdir_)
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
        ct_labvec_ = [];
        if nargin < 9
            C_indivcells_ = [];
            if nargin < 8
                k_ = 4;
                if nargin < 7
                    costfun_ = 'SumFit';
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

outstruct = struct;
sumfit_vec = zeros(1,length(ng_param_list_));
error_vec = sumfit_vec;

fprintf('Initializing preloaded gene indices\n')
tic
if length(ng_param_list_) > 5 && strcmp(method_,'MRx3') % heuristic criterion
    preloadinds = MRx3_Selector(genevct_,voxvgene_,max(ng_param_list_),lambda_);
else
    preloadinds = [];
end
toc

for i = 1:length(ng_param_list_)
    fprintf('Determining subset, nG parameter value %d/%d\n',i,length(ng_param_list_))
    param = ng_param_list_(i);
    % Create reduced versions of voxvgene and genevct according to the
    % method and parameter specified by the user
    tic
    if strcmp(costfun_,'SumFit')
        [E_red,C_red,nGen] = GeneSelector(genevct_,voxvgene_,gene_names_,param,lambda_,method_,preloadinds);
    else
        [E_red,C_red,nGen,~,C_ind_red] = GeneSelector_IndivCells(genevct_,voxvgene_,C_indivcells_,gene_names_,param,lambda_,method_,preloadinds);
    end
    toc
    % Infer cell density per voxel in arbitrary units
    fprintf('Nonnegative matrix inversion, nG parameter value %d/%d\n',i,length(ng_param_list_))
    tic
    B = CellDensityInference(E_red,C_red);
    toc
    outstruct(i).Bvals = B; 
    outstruct(i).nGen = nGen;
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

    % Calculate Pearson and Lin R
    [LinRstruct,PearsonStruct] = CorrelationsCalc(outstruct,i,matdir_);
    outstruct(i).LinR = LinRstruct;
    outstruct(i).Pearson = PearsonStruct;
    LinRnames = fieldnames(LinRstruct);
    Pnames = fieldnames(PearsonStruct);
    for j = 1:length(LinRnames)
        curparam_LinR(j,:) = LinRstruct.(LinRnames{j});
        curparam_Rval(j,:) = PearsonStruct.(Pnames{j});
    end
    LinR_pv(i,:) = curparam_LinR(1,:);
    LinR_sst(i,:) = curparam_LinR(2,:);
    LinR_vip(i,:) = curparam_LinR(3,:);
    LinR_all(i,:) = curparam_LinR(4,:);
    LinR_micro(i,:) = curparam_LinR(5,:);
    LinR_neuron(i,:) = curparam_LinR(6,:);
    Rval_pv(i,:) = curparam_Rval(1,:);
    Rval_sst(i,:) = curparam_Rval(2,:);
    Rval_vip(i,:) = curparam_Rval(3,:);
    Rval_all(i,:) = curparam_Rval(4,:);
    Rval_micro(i,:) = curparam_Rval(5,:);
    Rval_neuron(i,:) = curparam_Rval(6,:);

    % Calculate adjusted Kendall's tau for layer-type glutamatergic neurons
    fprintf('Calculating tau, nG parameter value %d/%d\n',i,length(ng_param_list_))
    ranks = [1 2 3 3 4 4 4];
    cell_inds = 9:15;
    cell_names = {'L2n3','L4','L5IT','L5PT','L6CT','L6IT','L6b'};
    taustruct = TauCalc(outstruct,i,cell_names,cell_inds,ranks,matdir_);
    outstruct(i).tau = taustruct.tau;

    % Calculate SumFit criterion
    fprintf('Calculating sum fit, nG parameter value %d/%d\n',i,length(ng_param_list_))
    sumfit = taustruct.tau + LinR_pv(i,3) + LinR_sst(i,3) + LinR_vip(i,3) + LinR_micro(i,3);
    sumfit_vec(i) = sumfit;
    outstruct(i).sumfit = sumfit;
    
    % Calculate classification error (SVM only)
    if strcmp(costfun_,'SVM')
        fprintf('Determining SVM classification error, nG parameter value %d/%d\n',i,length(ng_param_list_))
        tic
        ngroups = 1;
        options = statset('UseParallel',true); 
        errors = SVM_Classification_MISS(C_ind_red,ct_labvec_,k_,options,ngroups);
        outstruct(i).error = errors;
        toc
    elseif strcmp(costfun_,'LDA')
        fprintf('Determining LDA classification error, nG parameter value %d/%d\n',i,length(ng_param_list_))
        tic
        ngroups = k_;
        options = statset('UseParallel',true); 
        [errors,meanposts] = LDA_Classification_MISS(C_ind_red,ct_labvec_,options,ngroups);
        outstruct(i).error = errors;
        outstruct(i).posterior = meanposts;
        toc        
    end
    fprintf('Done, nG parameter value %d/%d\n',i,length(ng_param_list_))
end

fprintf('Determining optimal nG value\n');
if strcmp(costfun_,'SumFit')
    [~,optind] = max(sumfit_vec);
else
    nG_max = length(gene_names_);
    nG_rescaled  = (ng_param_list_ - min(ng_param_list_))/(nG_max - min(ng_param_list_));
    error_rescaled = (error_vec - min(error_vec))/max(error_vec - min(error_vec));
    sqdist = nG_rescaled.^2 + error_rescaled.^2;
    [~,optind] = min(sqdist);
end
end
