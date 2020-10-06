function voxvgene = ISH_Data_Extract_Zeisel(directory)

if nargin < 1
    directory = [cd filesep 'MatFiles'];
end
load([directory filesep 'ISH_data.mat'],'GENname','wherecoronal','GENsetid','GENsetid_coronal','V');
zeisel_genenames = h5read([directory filesep 'l5_all.agg.loom'],'/row_attrs/Gene');
zeisel_genenames = unique(sort(cellfun(@deblank, zeisel_genenames, 'UniformOutput', false)));

% Extract the gene names from the whole-brain expression data and create
% an index with which to sort the EntrezIDs in classstruct
uni_wb_genenames = unique(GENname(wherecoronal));
genebool = ismember(zeisel_genenames,uni_wb_genenames);
entrez_names = zeisel_genenames(genebool);

corinds = ismember(GENsetid,GENsetid_coronal);
genlist = GENname(corinds);
curnames = entrez_names;
gendex = ismember(genlist,curnames);
gendex = find(gendex);
percellmat = V(:,gendex);
redgenlist = genlist(gendex);
[redgenlist_sort,sortinds] = sort(redgenlist);
percellmat = percellmat(:,sortinds);
k = 1;
kt = length(redgenlist);
perCell = zeros(size(V,1), length(unique(redgenlist)));
for j = 1:length(unique(redgenlist_sort))
    if (k+1) <= kt && strcmp(redgenlist_sort{k},redgenlist_sort{k+1})
        if (k+2) <= kt && strcmp(redgenlist_sort{k+1},redgenlist_sort{k+2})
            if (k+3) <= kt && strcmp(redgenlist_sort{k+2},redgenlist_sort{k+3})
                perCell(:,j) = sum(percellmat(:,k:(k+3)),2)/4;
                k = k + 4;
            else
                perCell(:,j) = sum(percellmat(:,k:(k+2)),2)/3;
                k = k + 3;
            end
        else
            perCell(:,j) = sum(percellmat(:,k:(k+1)),2)/2;
            k = k + 2;   
        end
    else
        perCell(:,j) = percellmat(:,k);
        k = k + 1;
    end
end
voxvgene = perCell;
end
