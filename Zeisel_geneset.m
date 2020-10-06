metacols = h5read('/Users/christophermezias/Documents/MISS_General/MatFiles/l5_all.agg.loom','/col_attrs/ClusterName');
load('/Users/christophermezias/Documents/MISS_General/MatFiles/Zeisel_cellIDs.mat');
metacols = cellfun(@deblank,metacols,'UniformOutput',false);
[newcellnames,metainds] = sort(metacols);
cellIncl = logical(cellIncl);
inclcols = newcellnames(cellIncl);
exclcols = newcellnames(~cellIncl);
classkey = inclcols;

metacells = h5read('/Users/christophermezias/Documents/MISS_General/MatFiles/l5_all.agg.loom','/matrix');
cellcols = h5read('/Users/christophermezias/Documents/MISS_General/MatFiles/l5_all.loom','/col_attrs/ClusterName');
cellrows = h5read('/Users/christophermezias/Documents/MISS_General/MatFiles/l5_all.loom','/row_attrs/Gene');
newgennames = cellfun(@deblank,cellrows,'UniformOutput',false);
% newcellnames = cellfun(@deblank,cellcols,'UniformOutput',false);
[newgennames,geninds] = sort(newgennames);
metacells = metacells(metainds,geninds);
metacells = metacells(cellIncl,:);
[~,maxgen_percell] = max(metacells,[],2);
unimax = unique(maxgen_percell);
maxgenanmes = newgennames(unimax);

filedir = '/Users/christophermezias/Documents/MISS_General';
rawdata = readcell([filedir filesep 'mmc4.csv']);
cnames = rawdata(2:4:1058,1);
cnames{39} = 'OBDOP2';
genname_rows = 2:4:1058;
genname_raw = rawdata(genname_rows,3:end);
inclgen = ismember(cnames,classkey);
genname_incl = genname_raw(inclgen,:);
genname_vec = reshape(genname_incl,size(genname_incl,1)*size(genname_incl,2),1);
unigen_spex = unique(genname_vec);
genname_vec = [genname_vec;maxgenanmes];
unigen = unique(genname_vec);

load('/Users/christophermezias/Documents/MISS_General/MatFiles/ISH_gene_names.mat','ISH_gene_names');
entrez_inds = ismember(unigen,ISH_gene_names);
entrez_names = unigen(entrez_inds);

% reprows = [];
% for i = 1:length(entrez_names)
%     celltest = strcmp(genname_incl,entrez_names{i});
%     [rows,col] = find(celltest);
%     reprows = [reprows;rows];
% end
criteria = 1;
celltest = ismember(genname_incl,entrez_names);
gensum = sum(celltest,2);
reprows = (gensum>=criteria);
% reprows = unique(reprows);
newcnames = cnames(inclgen);
newcnames = newcnames(reprows);
repcellinds = ismember(classkey,newcnames);
repcells = classkey(repcellinds);
% check = strcmp(sort(newcnames),inclcols);
% exclrows = [15 71 165 167 171 172 191 196 198];
% newcellnames = cnames(inclgen);
% norepcells = newcellnames(exclrows);

Zeisel_gene_names = entrez_names;
save([filedir filesep 'MatFiles' filesep 'Zeisel_coronal_geneset.mat'],'Zeisel_gene_names','repcells','repcellinds');


