function [lambdastruct, peakind, sumfit_mat] = lambda_ParameterFitter(voxvgene, genevct, gene_names, ng_param_list, lambda_param_list)
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
if nargin < 5
    lambda_param_list = 50:50:500;
    if nargin < 4
        ng_param_list = 100:70:800;
        if nargin < 3
            error(sprintf('Error. \nUser must supply a voxel x gene matrix, a gene x cell type matrix, and a subsetting method'))
        end
    end
end 

sumfit_mat = zeros(length(ng_param_list),length(lambda_param_list));
for i = 1:length(lambda_param_list)
    fprintf('Running nG_ParameterFitter, lambda value %d/%d\n',i,length(lambda_param_list))
    lam = lambda_param_list(i);
    outstruct = nG_ParameterFitter(voxvgene,genevct,gene_names,'MRx3',ng_param_list,lam);
    if i == 1
        lambdastruct = outstruct;
    else
        lambdastruct = cat(2,lambdastruct,outstruct);
    end
    for j = 1:length(ng_param_list)
        sumfit_mat(j,i) = outstruct(j).sumfit;
    end
    fprintf('Done, lambda value %d/%d\n',i,length(lambda_param_list))
end
[~,peakind] = max(sumfit_mat(:));
end
