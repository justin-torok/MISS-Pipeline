metacols = h5read('/Users/christophermezias/Documents/l5_all.agg.loom','/col_attrs/ClusterName');
load('/Users/christophermezias/Documents/Zeisel_cellIDs.mat');
metacols = cellfun(@deblank,metacols,'UniformOutput',false);
[metacols,metainds] = sort(metacols);
cellIncl = logical(cellIncl);
inclcols = metacols(cellIncl);
exclcols = metacols(~cellIncl);
classkey = inclcols;

load('/Users/christophermezias/Documents/MISS_General/MatFiles/ISH_gene_names.mat','ISH_gene_names');

cellcols = h5read('/Users/christophermezias/Documents/l5_all.loom','/col_attrs/ClusterName');
cellrows = h5read('/Users/christophermezias/Documents/l5_all.loom','/row_attrs/Gene');
newgennames = cellfun(@deblank,cellrows,'UniformOutput',false);
newcellnames = cellfun(@deblank,cellcols,'UniformOutput',false);
[newgennames,geninds] = sort(newgennames);
[~,cellinds] = sort(newcellnames);
entrez_inds = ismember(newgennames,ISH_gene_names);
entrez_names = newgennames(entrez_inds);
gene_names = entrez_names;

metacells = h5read('/Users/christophermezias/Documents/l5_all.agg.loom','/matrix');
metacells = metacells(metainds,geninds);
newmetacells = metacells(cellIncl,entrez_inds);

groupinds = cell(length(metacols),1);
uniname = unique(newcellnames);
incl_len = 0;
for i = 1:length(uniname)
    curname = uniname{i};
    matchstrs = strcmp(newcellnames,curname);
    matchinds = find(matchstrs);
    groupinds(i) = {matchinds};
    if sum(strcmp(inclcols,curname))
        curlen = length(matchinds);
        incl_len = incl_len + curlen;
    end
end

C_indivcells = zeros(incl_len,length(entrez_names));
ct_labvec = zeros(incl_len,1);
ct_namevec = cell(incl_len,1);
dex = 0;
id_dex = 1;
for j = 1:length(groupinds)
    curname = uniname{j};
    if sum(strcmp(inclcols,curname))
        curinds = groupinds{j};
        read_start = [curinds(1) 1];
        read_count = [length(curinds) length(newgennames)];
        hdf_indcells = h5read('/Users/christophermezias/Documents/l5_all.loom','/matrix',read_start,read_count);
%         if j == 10
%             genname_test = h5read('/Users/christophermezias/Documents/l5_all.loom','/row_attrs/Gene');
%             genname_test = cellfun(@deblank,genname_test,'UniformOutput',false);
%             genname_test = genname_test(geninds);
%             genname_test = genname_test(entrez_inds);
%         end
        hdf_indcells = hdf_indcells(:,geninds);
        hdf_indcells = hdf_indcells(:,entrez_inds);
        C_indivcells(dex+1:dex+length(curinds),:) = hdf_indcells;
        ct_labvec(dex+1:dex+length(curinds),1) = ones(length(curinds),1) * id_dex;
        ct_namevec(dex+1:dex+length(curinds),1) = newcellnames(curinds(1));
        dex = dex + length(curinds);
        id_dex = id_dex + 1;
        clear hdf_indcells
    end
end

dup_inds = zeros(1,length(gene_names));
for k = 1:length(gene_names)
    curname = gene_names{k};
    isduplicate = strcmp(gene_names,curname);
    isduplicate = find(isduplicate);
    if length(isduplicate) > 1
        curgens = newmetacells(:,isduplicate);
        mnexp = mean(curgens,2);
        newmetacells(:,isduplicate(1)) = mnexp;
        curgens_indc = C_indivcells(:,isduplicate);
        mn_indcells = mean(curgens_indc,2);
        C_indivcells(:,isduplicate(1)) = mn_indcells;
        if length(isduplicate) > 2
            dup_inds(isduplicate(2:end)) = 1;
        elseif length(isduplicate) == 2
            dup_inds(isduplicate(2)) = 1;
        end
    end
end
dup_inds = logical(dup_inds);
meanexprmat = newmetacells;
meanexprmat(:,dup_inds) = [];
entrez_names(dup_inds) = [];
C_indivcells(:,dup_inds) = [];

clearvars -except C_indivcells ct_labvec ct_namevec entrez_names meanexprmat classkey
        
save('/Users/christophermezias/Documents/MISS_General/MatFiles/Zeisel_extract.mat','-v7.3');
    

% facs = factor(length(cellcols));
% facs = [facs(end) length(cellcols)/facs];
% dat_divisor = max(facs);
% nreps = min(facs);
% 
% dex = 0;
% C_indcells = zeros(length(cellcols),length(entrez_names));
% for i = 1:nreps
%     read_start = [dex+1 1];
%     read_count = [dat_divisor length(cellrows)];
%     hdf_indcells = h5read('/Users/christophermezias/Documents/l5_all.loom','/matrix',read_start,read_count);
%     hdf_indcells = hdf_indcells(cellinds(dex+1:dex+dat_divisor),geninds);
%     hdf_indcells = hdf_indcells(:,entrez_inds);
%     
%     C_indcells(dex+1:dex+dat_divisor,:) = hdf_indcells;
%     dex = dex + dat_divisor;
%     
%     
%     
%     
%     
% 
% 
% 
% 
% 
