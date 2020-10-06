function [meanexprmat, classkey, entrez_names] = Cell_Type_Data_Extract_Tasic(classstruct,directory)

if nargin < 2
    directory = [cd filesep 'MatFiles'];
end

load([directory filesep 'ISH_data.mat'],'GENname','wherecoronal');
load([directory filesep 'mouse_VISp_2018-06-14_genes-rows.mat'],'filestruct'); 
genesrows = filestruct;

% Extract the gene names from the whole-brain expression data and create
% an index with which to sort the EntrezIDs in classstruct
uni_wb_genenames = unique(GENname(wherecoronal));
genebool = ismember(genesrows.gene_symbol,uni_wb_genenames);
ct_entrez = genesrows.gene_entrez_id(genebool);
entrezbool = ismember(classstruct.VISp.EntrezID,ct_entrez);
entrez_names = genesrows.gene_symbol(entrezbool);

regions = {'VISp','ALM','LGd'};
classkey = cell(2,1);
exprmat = [];
cellcountvec = [];
endoexpr = []; astrexpr = []; oligexpr = []; micrexpr = [];
endocount = []; astrcount = []; oligcount = []; micrcount = [];
for i=1:length(regions)
    classfields = fieldnames(classstruct.(regions{i}));
    classfields = setdiff(classfields,{'EntrezID','Low_Quality'});
    for j=1:length(classfields)
        subclassfields = fieldnames(classstruct.(regions{i}).(classfields{j}));
        subclassfields = setdiff(subclassfields,{'truename','meanexpr','cellcount'});
        if ismember(classfields{j},{'GABAergic', 'Glutamatergic'})
            for k=1:length(subclassfields)
                meanexpr = classstruct.(regions{i}).(classfields{j}).(subclassfields{k}).meanexpr(entrezbool);
                cellcount = classstruct.(regions{i}).(classfields{j}).(subclassfields{k}).cellcount;
                expr = meanexpr * cellcount;
                if size(classkey,2) == 1
                    duplicatesubclass = 0;
                else
                    duplicatesubclass = zeros(1,(size(classkey,2)-1));
                    for l = 1:length(duplicatesubclass)
                        subclasstest = classkey{2,l+1};
                        duplicatesubclass(l) = strcmp(subclassfields{k},subclasstest);
                    end
                end
                
                if sum(duplicatesubclass)==0
                    classkey(:,end+1) = cell(2,1);
                    classkey{1,end} = classfields{j};
                    classkey{2,end} = subclassfields{k};
                    if isempty(exprmat)
                        exprmat = expr;
                        cellcountvec = cellcount;
                    else
                        exprmat = cat(2,exprmat,expr);
                        cellcountvec = [cellcountvec, cellcount];
                    end                    
                else
                    dupind = logical(duplicatesubclass);
                    exprmat(:,dupind) = exprmat(:,dupind) + expr;
                    cellcountvec(dupind) = cellcountvec(dupind) + cellcount;
                end
            end
        elseif ismember(classfields{j},{'Non_Neuronal'})
            for k=1:length(subclassfields)
                meanexpr = classstruct.(regions{i}).(classfields{j}).(subclassfields{k}).meanexpr(entrezbool);
                cellcount = classstruct.(regions{i}).(classfields{j}).(subclassfields{k}).cellcount;
                expr = meanexpr * cellcount;
                if strcmp(subclassfields{k}(1:3),'Ast')
                    if isempty(astrexpr)
                        astrexpr = expr;
                        astrcount = cellcount;
                    else
                        astrexpr = astrexpr + expr;
                        astrcount = astrcount + cellcount;
                    end
                elseif strcmp(subclassfields{k}(1:3),'Oli')
                    if isempty(oligexpr)
                        oligexpr = expr;
                        oligcount = cellcount;
                    else
                        oligexpr = oligexpr + expr;
                        oligcount = oligcount + cellcount;
                    end            
                elseif strcmp(subclassfields{k}(1:3),'Mac')
                    if isempty(micrexpr)
                        micrexpr = expr;
                        micrcount = cellcount;
                    else
                        micrexpr = micrexpr + expr;
                        micrcount = micrcount + cellcount;
                    end
                end
            end
        elseif ismember(classfields{j},{'Endothelial'})
            meanexpr = classstruct.(regions{i}).(classfields{j}).meanexpr(entrezbool);
            cellcount = classstruct.(regions{i}).(classfields{j}).cellcount;
            expr = meanexpr * cellcount;
            if isempty(endoexpr)
                endoexpr = expr;
                endocount = cellcount;
            else
                endoexpr = endoexpr + expr;
                endocount = endocount + cellcount;
            end
        end
    end
end
exprmat = cat(2, exprmat, astrexpr, micrexpr, oligexpr, endoexpr);
cellcountvec = cat(2, cellcountvec, astrcount, micrcount, oligcount, endocount);
meanexprmat = zeros(size(exprmat));
for i = 1:size(exprmat,2)
    meanexprmat(:,i) = exprmat(:,i) / cellcountvec(i);
end

classkey = classkey(2,2:end);
classkey = cat(2,classkey,{'Astro','Micro','Oligo','Endo'});

end