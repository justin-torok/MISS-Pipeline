function [lambda_opt,sigma_opt,ngstar_opt] = Hyperparameter_Optimizer(method,...
                                                nG_vals, lambda_list,...
                                                sigma_list, k, cv, directory)

% This function outputs the lambda and sigma hyperparameter combination
% that produces an nG* value that both minimizes the nearest centroid
% clustering-based loss function and maximizes the fit to expectations via
% sumfit, which it determines by performing a grid search over supplied
% ranges of each set of values.

if nargin < 7
    directory = [cd filesep 'MatFiles'];
    if nargin < 6
        cv = 1;
        if nargin < 5
            k = 5;
            if nargin < 4
                sigma_list = 3000:25:6000;
                if nargin < 3
                    lambda_list = 0:50:350;
                    if nargin < 2
                        if nargin < 1
                            error(sprintf('Error. \nUser must supply a subsetting method'))
                        end
                        if ismember(method, {'MRx3','mRMR','colAMD'})
                            nG_vals = [50:50:300,301:1200,1255:50:2155,2255:100:3655];
                        elseif strcmp(method, 'DBSCAN')
                            nG_vals = [0.0001, 0.0004:0.0001:0.001, 0.002:0.001:0.005, 0.015:0.01:0.195]; 
                        elseif strcmp(method, 'Entropy')
                            nG_vals = [0.25,0.5:0.125:2.75,2.8:0.05:3.25];
                        elseif strcmp(method, 'Zeisel')
                            nG_vals = [2:5, 6:2:50];
                        else
                            error('Error. \n%s is an incorrect subsetting method identifier',method)
                        end
                    end    
                end
            end
        end
    end
end

load([directory filesep 'Tasic_Inputs.mat'],'voxvgene','genevct',...
    'gene_names','ct_labvec','C_indivcells');

sumfits = zeros(length(lambda_list),length(sigma_list)); ngstars = sumfits;
for i = 1:length(lambda_list)
    fprintf('Determining sumfit at nG*, lambda parameter value %d/%d\n',i,length(lambda_list))
    lambda = lambda_list(i);
    fitstruct = nG_ParameterFitter(voxvgene, genevct, gene_names,...
                                   method, C_indivcells, ct_labvec,...
                                   nG_vals, lambda, 4400, k, cv, directory);    
    [sumfitvec,minngvec] = SumFit_Optimizer(fitstruct,method,sigma_list,directory);
    sumfits(i,:) = sumfitvec; ngstars(i,:) = minngvec;
    fprintf('Done, lambda parameter value %d/%d\n',i,length(lambda_list))
end
lambdamat = repmat(lambda_list.',1,length(sigma_list));
sigmamat = repmat(sigma_list,length(lambda_list),1);
[~,maxind] = max(sumfits(:));
lambda_opt = lambdamat(maxind); sigma_opt = sigmamat(maxind);
ngstars = ngstars(:); ngstar_opt = ngstars(maxind);
end