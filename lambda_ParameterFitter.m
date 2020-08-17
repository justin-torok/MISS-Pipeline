function [lambdastruct, peakind] = lambda_ParameterFitter(voxvgene, genevct, gene_names, ng_param_list, lambda_param_list, costfun, k, C_indivcells, ct_labvec, matdir)
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
% different aspects of the data processing and quantitative assessment. 

% Default lambda and nG parameter lists
if nargin < 10
    matdir = cd;
    if nargin < 9
        ct_labvec = [];
        if nargin < 8
            C_indivcells = [];
            if nargin < 7
                k = 4;
                if nargin < 6
                    costfun = 'SumFit';
                    if nargin < 5
                        lambda_param_list = 50:50:500;
                        if nargin < 4
                            ng_param_list = 100:70:800;
                            if nargin < 3
                                error(sprintf('Error. \nUser must supply a voxel x gene matrix, a gene x cell type matrix, and a subsetting method'))
                            end
                        end
                    end 
                end
            end
        end
    end
end

cost_mat = zeros(length(ng_param_list),length(lambda_param_list));
for i = 1:length(lambda_param_list)
    fprintf('Running nG_ParameterFitter, lambda value %d/%d\n',i,length(lambda_param_list))
    lam = lambda_param_list(i);
    outstruct = nG_ParameterFitter(voxvgene,genevct,gene_names,'MRx3',ng_param_list,lam,costfun,k,C_indivcells,ct_labvec,matdir);
    if i == 1
        lambdastruct = outstruct;
    else
        lambdastruct = cat(2,lambdastruct,outstruct);
    end
    for j = 1:length(ng_param_list)
        if strcmp(costfun,'SumFit')
            cost_mat(j,i) = -outstruct(j).sumfit;
        else
            cost_mat(j,i) = mean(outstruct(j).error);
        end
    end
    fprintf('Done, lambda value %d/%d\n',i,length(lambda_param_list))
end
[~,peakind] = min(cost_mat(:));
end
