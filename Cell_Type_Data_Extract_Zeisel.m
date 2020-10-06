function [genevct, C_indivcells, classkey, entrez_names, ct_labvec, ct_namevec] = Cell_Type_Data_Extract_Zeisel(directory)

if nargin < 1
    directory = [cd filespe 'MatFiles'];
end

load([directory filesep 'Zeisel_cellIDs.mat'],'cellIncl');
load([directory filesep 'ISH_gene_names.mat'],'ISH_gene_names');
metacols = h5read([directory filesep 'l5_all.agg.loom'],'/col_attrs/ClusterName');
metacols = cellfun(@deblank,metacols,'UniformOutput',false);
[metacols,metainds] = sort(metacols);
cellIncl = logical(cellIncl);
inclcols = metacols(cellIncl);
% exclcols = metacols(~cellIncl);
classkey = inclcols;

cellcols = h5read([directory filesep 'l5_all.loom'],'/col_attrs/ClusterName');
cellrows = h5read([directory filesep 'l5_all.loom'],'/row_attrs/Gene');
newgennames = cellfun(@deblank,cellrows,'UniformOutput',false);
newcellnames = cellfun(@deblank,cellcols,'UniformOutput',false);
[newgennames,geninds] = sort(newgennames);
% [~,cellinds] = sort(newcellnames);
entrez_inds = ismember(newgennames,ISH_gene_names);
entrez_names = newgennames(entrez_inds);
gene_names = entrez_names;

metacells = h5read([directory filesep 'l5_all.agg.loom'],'/matrix');
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
        hdf_indcells = h5read([directory filesep 'l5_all.loom'],'/matrix',read_start,read_count);
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
genevct = newmetacells;
genevct(:,dup_inds) = [];
entrez_names(dup_inds) = [];
C_indivcells(:,dup_inds) = [];

genevct = genevct.'; ct_labvec = ct_labvec.'; ct_namevec = ct_namevec.';
C_indivcells = C_indivcells.'; classkey = classkey.';

end