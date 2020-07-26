function [E_red,C_red,nGen,reduced_gene_names,C_ind_red] = GeneSelector_IndivCells(genevct,voxvgene,C_indivcells,gene_names,ngen_param,lambda,method,preloadinds)
% Function that performs various types of subset selection implemented in
% Mezias et al, 2020, including the ideal method MRx3. The user MUST supply
% the following inputs:
% genevct - a scRNAseq-based expression matrix that is genes x cell types
% in dimension (output of ProcessedData_Generator.m)
% voxvgene - an ISH-based expression matrix that is voxels x genes in
% dimension (output of ProcessedData_Generator.m)
% gene_names - a genes x 1 cell array of gene names (output of 
% ProcessedData_Generator.m) 
% ngen_param - integer, definition depends on subsetting method For mRMR,
% MRx3, and colAMD, ngen_param is the size of the target gene set. For 
% Entropy, it is  the threshold entropy value above which genes are not
% selected. For DBSCAN, it is the value of the neighborhood search radius. 
% For Zeisel, it is the number of genes to select per cell type. 
%
% The user may also choose to supply:
% lambda - MRx3-specific parameter that controls the relative weighting of
% gene projection error; if not provided, set to 150
% method - string specifying the subset selection method; if not provided,
% set to MRx3
% preloadinds - a 1 x length(gene_names) numeric array of indices for the
% gene order of a given subset method. It is more computationally efficient
% for parameter sweeps of MRx3 to determine the order of the genes selected
% across all genes in gene_names and then subset based on the first
% ngen_param entries of the list, as opposed to calling MRx3_Selector many
% times. This is unnecessary but would still work in principle for any of
% the other gene selection methods as well. In this case the ngen_param the
% user supplies should serve the same function as if they were performing
% MRx3, mRMR, or colAMD.
%
% The outputs are:
% E_red - reduced genes x voxels expression matrix
% C_red - reduced genes x cell types expression matrix
% nGen - number of genes in the reduced set
% reduced_gene_names - reduced genes x 1 cell array of gene names for the 
% reduced set

% Default option is MRx3 with lambda = 150

if nargin < 8
    preloadinds = [];
    if nargin < 7
        method = 'MRx3';
        if nargin < 6
            lambda = 150;
        end
    end
end


if isempty(preloadinds)
    % Normalizations for genevct that are necessary for performing several subset
    % selection methods (MRx3 does these internally)
    ctmean = mean(genevct,1);
    ctmean = repmat(ctmean,size(genevct,1),1);
    ctnorm = genevct ./ ctmean;
    ctnorm(ctnorm<(0.1*std(nonzeros(ctnorm)))) = 0;

    genesum = sum(ctnorm,2);
    genesum = repmat(genesum,1,size(genevct,2));
    genenorm = ctnorm ./ genesum;
    genenorm(isnan(genenorm)) = 0;

    % Subset selection performed according to the "method" argument above. If
    % not supplied, it is assumed to be MRx3
    if strcmp(method,'MRx3')
        mrmr_n = ngen_param;
        mrmrinds = MRx3_Selector(genevct,voxvgene,mrmr_n,lambda);
        mrmrinds = sort(mrmrinds);
        reduced_gene_names = gene_names(mrmrinds);
    elseif strcmp(method,'mRMR')
       mrmr_n = ngen_param;
       mrmrinds = mRMR_Selector(ctnorm,mrmr_n,'Quo');
       reduced_gene_names = gene_names(mrmrinds);
       reduced_gene_names = unique(reduced_gene_names);
    elseif strcmp(method,'DBSCAN')
        reduced_gene_names = 'a';
        epsilon = ngen_param; % distance parameter for DBSCAN
        for i = 1:size(genenorm,2)
            curtype = genenorm(:,i);
            dbinds = dbscan(curtype,epsilon,30);
            genenames{i} = gene_names(dbinds == -1);
            genenames{i} = sort(genenames{i});
            reduced_gene_names = [reduced_gene_names;genenames{i}];        
        end
        reduced_gene_names(1) = [];
        reduced_gene_names = unique(reduced_gene_names);
    elseif strcmp(method,'Entropy')
        entthresh = ngen_param; % entropy score threshold
        entropyscore = zeros(1,size(genenorm,1));
        for i = 1:size(genenorm,1)
            curtype = genenorm(i,:);
            curtype(curtype == 0) = 1; % necessary for entropy calculation
            entropyscore(i) = -dot(log(curtype),curtype);
        end
        entropyscore(entropyscore==0) = max(entropyscore);
        entinds = find(entropyscore < entthresh);
        reduced_gene_names = gene_names(entinds);
    elseif strcmp(method,'colAMD')
        thresh = 0.04; %making matrix x% sparse before colamd
        normthresh = genenorm;
        normthresh(normthresh<thresh) = 0;
        normthresh = normthresh.';
        CC = colamd(normthresh);
        cainds = CC(1:ngen_param);
        reduced_gene_names = gene_names(cainds);
    elseif strcmp(method,'Zeisel')
        reduced_gene_names = 'a';
        for i = 1:size(genenorm,2)
            [~, maxind] = max(ctnorm(:,i));
            maxname = gene_names(maxind);
            [~,inds_new] = sort(genenorm(:,i),'descend');
            noabovethresh = ngen_param-1;
            specnames = gene_names(inds_new(1:noabovethresh));
            genenames{i} = [maxname; specnames];
            genenames{i} = sort(genenames{i});
            reduced_gene_names = [reduced_gene_names;genenames{i}];
        end
        reduced_gene_names(1) = [];
        reduced_gene_names = unique(reduced_gene_names);
    end
else % User has supplied preloadinds
    redinds = sort(preloadinds(1:ngen_param));
    reduced_gene_names = gene_names(redinds);  
end
% generating indices for diff exp genes, choosing in spatexp and genevct
gendex = ismember(gene_names,reduced_gene_names);
gendex = find(gendex);
E_red = (voxvgene(:,gendex)).';
C_red = genevct(gendex,:);
C_ind_red = C_indivcells(gendex,:);
nGen = size(C_red,1);

end