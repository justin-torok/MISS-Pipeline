function [C_indivcells, classkey, entrez_names, ct_labvec] = Cell_Type_Data_Extract_IndivCells_Tasic(classstruct,directory)

if nargin < 2
    directory = [cd filesep 'MatFiles'];
end

load([directory filesep 'ISH_data.mat'],'GENname','wherecoronal');
load([directory filesep 'mouse_VISp_2018-06-14_genes-rows.mat'],'filestruct'); 
genesrows = filestruct;

uni_wb_genenames = unique(GENname(wherecoronal));
genebool = ismember(genesrows.gene_symbol,uni_wb_genenames);
ct_entrez = genesrows.gene_entrez_id(genebool);
entrezbool = ismember(classstruct.VISp.EntrezID,ct_entrez);
entrez_names = genesrows.gene_symbol(entrezbool);
members = ismember(filestruct.gene_symbol,entrez_names);

C_indivcells = [];
regnames = fieldnames(classstruct);
sublabs = {};
typelabs = {};

for i = 1:length(regnames)
    classnames = fieldnames(classstruct.(regnames{i}));
    exclinds = ismember(classnames,{'EntrezID','Outlier','Low_Quality'});
    classnames(exclinds) = [];
    for j = 1:length(classnames)
        typenames = fieldnames(classstruct.(regnames{i}).(classnames{j}));
        typenames([1 end-1:end]) = [];
        for k = 1:length(typenames)
            curtype = typenames{k};
            subnames = fieldnames(classstruct.(regnames{i}).(classnames{j}).(typenames{k}));
            subnames([1 end-1:end]) = [];
            for h = 1:length(subnames)
                C_indivcells = [C_indivcells classstruct.(regnames{i}).(classnames{j}).(typenames{k}).(subnames{h}).normalizedexpr(members,:)];
                cursublabs = cell(1,classstruct.(regnames{i}).(classnames{j}).(typenames{k}).(subnames{h}).cellcount);
                cursub = subnames{h};
                cursublabs(:) = {cursub};
                sublabs = [sublabs cursublabs];
                curtypelabs = cell(size(cursublabs));
                curtypelabs(:) = {curtype};
                typelabs = [typelabs curtypelabs];
            end
        end
    end
end

for i = 1:length(typelabs)
    curtype = typelabs{i};
    if strcmp(curtype,'Olig1')
        typelabs{i} = 'Oligo';
    elseif strcmp(curtype,'OPC')
        typelabs{i} = 'Oligo';
    elseif strcmp(curtype,'Macrophage')
        typelabs{i} = 'Macro';
    elseif strcmp(curtype,'Micro')
        typelabs{i} = 'Macro';
    elseif strcmp(curtype,'Peri')
        typelabs{i} = 'Endo';
    elseif strcmp(curtype,'SMC')
        typelabs{i} = 'Endo';
    end
end

load([directory filesep 'PresetInputs.mat'],'classkey');
ct_labvec = zeros(size(typelabs));
subt_labvec = zeros(size(sublabs));
dex = 1;
for i = 1:length(classkey)
    membs = ismember(typelabs,classkey{i});
    ct_labvec(membs) = i;
    cursubkeys = unique(sublabs(membs));
    for j = 1:length(cursubkeys)
        submembs = ismember(sublabs,cursubkeys{j});
        subt_labvec(submembs) = dex;
        dex = dex + 1;
    end 
end

excl_type = 0;

excltype_inds = (ct_labvec==excl_type);
C_indivcells(:,excltype_inds) = [];
ct_labvec(excltype_inds) = [];
subt_labvec(excltype_inds) = [];
typelabs(excltype_inds) = [];
if sum(excl_type) > 0
    classkey(excl_type) = [];
end

clearvars -except C_indivcells entrez_names classkey ct_labvec

end

