function [outstruct, peakind] = nG_ParameterFitter(voxvgene, genevct, gene_names, method, ng_param_list, lambda, costfun, k, C_indivcells, ct_labvec, matdir)
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
    matdir = cd;
    if nargin < 10
        ct_labvec = [];
        if nargin < 9
            C_indivcells = [];
            if nargin < 8
                k = 4;
                if nargin < 7
                    costfun = 'SumFit';
                    if nargin < 6
                        lambda = 150;
                        if nargin < 5
                            if nargin < 4
                                error(sprintf('Error. \nUser must supply a voxel x gene matrix, a gene x cell type matrix, a gene names cell array, and a subsetting method'))
                            end
                            if ismember(method, {'MRx3','mRMR','colAMD'})
                                ng_param_list = [100:70:1500, 1750:250:3500, 3855]; 
                            elseif strcmp(method, 'DBSCAN')
                                ng_param_list = [0.0001, 0.0004:0.0001:0.001, 0.002:0.001:0.005, 0.015:0.01:0.195]; 
                            elseif strcmp(method, 'Entropy')
                                ng_param_list = [0.25,0.5:0.125:2.75,2.8:0.05:3.25];
                            elseif strcmp(method, 'Zeisel')
                                ng_param_list = [2:5, 6:2:50]; % made up something here
                            else
                                error('Error. \n%s is an incorrect subsetting method identifier',method)
                            end
                        end    
                    end
                end
            end
        end
    end
end
outstruct = struct;
sumfit_vec = zeros(1,length(ng_param_list));
fprintf('Initializing preloaded gene indices\n')
if (length(ng_param_list) > 5 || strcmp(costfun,'SVM')) && strcmp(method,'MRx3') % heuristic criterion
    preloadinds = MRx3_Selector(genevct,voxvgene,length(gene_names),lambda);
else
    preloadinds = [];
end

for i = 1:length(ng_param_list)
    fprintf('Determining subset, nG parameter value %d/%d\n',i,length(ng_param_list))
    param = ng_param_list(i);
    % Create reduced versions of voxvgene and genevct according to the
    % method and parameter specified by the user
    if strcmp(costfun,'SumFit')
        [E_red,C_red,nGen] = GeneSelector(genevct,voxvgene,gene_names,param,lambda,method,preloadinds);
    else
        [E_red,C_red,nGen,~,C_ind_red] = GeneSelector_IndivCells(genevct,voxvgene,C_indivcells,gene_names,param,lambda,method,preloadinds);
    end
    
    % Infer cell density per voxel in arbitrary units
    fprintf('Nonnegative matrix inversion, nG parameter value %d/%d\n',i,length(ng_param_list))
    B = CellDensityInference(E_red,C_red);
    outstruct(i).Bvals = B; 
    outstruct(i).nGen = nGen;
    if strcmp(method,'MRx3')
        outstruct(i).lambda = lambda;
    end

    % Convert arbitrary densities to counts per voxel
    Bcorrected = Density_to_Counts(B,matdir);
    outstruct(i).corrB = Bcorrected;

    % Sum and average over CCF regions
    [sumB,meanB] = Voxel_To_Region(Bcorrected,matdir);
    outstruct(i).Bsums = sumB; % total cells per region
    outstruct(i).Bmeans = meanB; % mean cell count per region

    % Calculate Pearson and Lin R
    [LinRstruct,PearsonStruct] = CorrelationsCalc(outstruct,i,matdir);
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
    fprintf('Calculating tau, nG parameter value %d/%d\n',i,length(ng_param_list))
    ranks = [1 2 3 3 4 4 4];
    cell_inds = 9:15;
    cell_names = {'L2n3','L4','L5IT','L5PT','L6CT','L6IT','L6b'};
    taustruct = TauCalc(outstruct,i,cell_names,cell_inds,ranks,matdir);

    % Calculate SumFit criterion
    fprintf('Calculating sum fit, nG parameter value %d/%d\n',i,length(ng_param_list))
    sumfit = taustruct.tau + LinR_pv(i,3) + LinR_sst(i,3) + LinR_vip(i,3) + LinR_micro(i,3);
    sumfit_vec(i) = sumfit;
    outstruct(i).sumfit = sumfit;
    
    % Calculate classification error (SVM only)
    if strcmp(costfun,'SVM')
        fprintf('Determining SVM classification error, nG parameter value %d/%d\n',i,length(ng_param_list))
        x = (C_ind_red - mean(C_ind_red, 1) ./ std(C_ind_red, 1))';
        y = ct_labvec.';
        options = statset('UseParallel',true);
        Mdl = fitcecoc(x,y,'Options',options);
        rng(1); % set seed
        CVMdl = crossval(Mdl, 'Options', options, 'Kfold', k);
        error_vec(i) = kfoldLoss(CVMdl, 'Options', options);
        outstruct(i).error = error_vec(i);
    end
    fprintf('Done, nG parameter value %d/%d\n',i,length(ng_param_list))
end

if strcmp(costfun,'SumFit')
    [~,peakind] = max(sumfit_vec);
else
    [~,peakind] = min(error_vec);
end
end
