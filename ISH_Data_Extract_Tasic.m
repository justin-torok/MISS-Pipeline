function regvgene = ISH_Data_Extract_Tasic(classstruct, directory)

if nargin < 2
    directory = [cd filesep 'MatFiles'];
end
load([directory filesep 'ISH_data.mat'],'GENname','wherecoronal','GENsetid','GENsetid_coronal','V');
load([directory filesep 'mouse_VISp_2018-06-14_genes-rows.mat'],'filestruct'); 
genesrows = filestruct;

% Extract the gene names from the whole-brain expression data and create
% an index with which to sort the EntrezIDs in classstruct
uni_wb_genenames = unique(GENname(wherecoronal));
genebool = ismember(genesrows.gene_symbol,uni_wb_genenames);
ct_entrez = genesrows.gene_entrez_id(genebool);
entrezbool = ismember(classstruct.VISp.EntrezID,ct_entrez);
entrez_names = genesrows.gene_symbol(entrezbool);

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
regvgene = perCell;
end
