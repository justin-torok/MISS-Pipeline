function [genevct,voxvgene,ctkey,icell_key,genes] = ProcessedData_Generator(user_indivcells,user_labvec,user_gennames,userkey,user_weightvec,fpath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constructing appropriate dataset for Linear Inference Cell-Type Mapping

% Preset MatFile Input Variables: 
% meanexprmat = matrix of weighted mean scRNAseq data, with weights 
% assigned by correlation with the centroid of each cell-type cluster.

% ct_labvec = a numeric label vector from 1:25 for current cell types that 
% is n-single cells long. 

% classkey = string array of cell-type names, 25 elements long.

% regvgene = a voxel X gene matrix of ISH expression in ABI ISH mouse atlas 
% space.

% entrez_names = a string vector of entrez gene names in common between the 
% ISH and scRNAseq datasets.

% User Defined Input Variables: 
% user_indivcells = a matrix of cell-types the user supplies for mapping 
% additional cell-types, genes X n-cells.

% user_labvec: a numeric label vector that is 1 X n-cells from user data 
% long, akin to ct_labvec for the preset dataset. If a cell is from a type 
% overlapping with one in the preset dataset, then that cell-type position
% in the string vector *classkey* is it's numeric ID. All cells not from a 
% type found in the preset dataset or classkey should be given labels of 
% 26:26+n-cells according to the larger type-group each individual cell is 
% drawn from.

% userkey: a string array of cell type names unique to the user data,
% input in the same order as numeric labels 26:26+n-cells from user
% defined data.

% user_gennames: a string vector of entrez gene names found in the
% user's dataset.

% user_weightvec: a 1 X n-cells or n-cells X 1 vector of weights for each
% cell assigned by that individual cells correlation with the centroid of
% it's cell-type group/cluster, akin to the weighting done to meanexprmat
% from the preset data above. If not specified by the user, this is
% assigned to be a vector of 1s.

% fpath: the filepath telling the computer where to draw preset data from.
% If this is not assigned by the user, this is set to the CD.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Defining filepath as cd if not user input
if nargin < 6
    fpath = cd;
end

%loads C_indivcells, ct_labvec, classkey, regvgene, entrez_names
load([fpath filesep 'PresetInputs.mat']); 
entrez_names = sort(entrez_names);

if nargin >= 1
    
    %defining weighting vector for CT averaging if not specified by user
    if nargin < 5
        user_weightvec = ones(1,size(user_indivcells,2));
    end
    
    %sorting and getting unity between preset and user defined genes
    user_gennames = sort(user_gennames);
    unityinds_user = ismember(user_gennames,entrez_names);
    unitynames = user_gennames(unityinds_user);
    unityinds_pre = ismember(entrez_names,user_gennames);
    areoverlap = ismember(userkey,classkey);
    sumoverlap = sum(areoverlap);
    uniquekey = userkey;
    meancells = zeros(length(unitynames),length(classkey)+length(userkey)-sumoverlap);
    uniquekey(areoverlap) = [];
    
    %averaging per-cell type among types in preset data, including user
    %data among preset cell-types
    for i = 1:length(classkey)
        curct_pre = classkey{i};
        whereoverlap = find(ismember(userkey,curct_pre));
        cur_userinds = (user_labvec==whereoverlap);
        cur_precells = meanexprmat(unityinds_pre,i);
        if ~isempty(cur_userinds)
            cur_usercells = user_indivcells(unityinds_user,cur_userinds);
            cur_weights = user_weightvec(cur_userinds);
            ncells_user = length(cur_weights);
            if size(cur_weights,1) == 1
                cur_weights = cur_weights.';
            end
            curmn_usercells = (cur_usercells*cur_weights) / sum(cur_weights);
            ncells_pre = length(find(ct_labvec==i));
            curmncells = (cur_precells*ncells_pre + curmn_usercells*ncells_user) / (ncells_pre + ncells_user);
        else
            curmncells = cur_precells;
        end
        meancells(:,i) = curmncells;
    end
    
    %averaging per cell-type data for user supplied cells only
    if isempty(uniquekey)
        newkey = classkey;
    else
        uniquemninds = (length(classkey)+1):(length(classkey)+length(uniquekey));
        whereunique = find(ismember(userkey,uniquekey));
        for j = 1:length(uniquekey)
            cur_userinds = (user_labvec==whereunique(j));
            cur_usercells = user_indivcells(unityinds_user,cur_userinds);
            cur_weights = user_weightvec(cur_userinds);
            if size(cur_weights,1)==1
                cur_weights = cur_weights.';
            end
            curmncells = (cur_usercells*cur_weights)/sum(cur_weights);
            meancells(:,uniquemninds(j)) = curmncells;
        end
        if size(uniquekey,1) > 1
            uniquekey = uniquekey.';
        end
        newkey = [classkey uniquekey];
    end
    
    newlabvec = [ct_labvec user_labvec];
    unity_ish = regvgene(:,unityinds_pre);
    
else
    
    %averaging per cell-type data if only using preset data
    meancells = meanexprmat;
    newlabvec = ct_labvec;
    newkey = classkey;
    unitynames = entrez_names;
    unity_ish = regvgene;
    
end

%decalring output variables
voxvgene = unity_ish;
ctkey = newkey;
icell_key = newlabvec;
genes = unitynames;
genevct = meancells;

end



    




