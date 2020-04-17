function metric_struct = TauCalc(outstruct,idx,cell_names,cell_inds,ranks,directory)
% This function calculates the adjusted Kendall's tau correlation for
% layer-specific glutamatergic neurons in the neocortex, which is one of
% the components of sum fit. The algorithm determines how well the seven
% neuronal types are segregated in their expected layers by determining the
% midline of their cell density banding patterns and then ranking the
% distances of these midlines from the cortical surface. The tau value is
% calculated across all cortical slices containing sufficient neocortex and
% is penalized by neuronal types with insufficient neocortical signal
% relative to background.
%
% The user must supply an outstruct and a paramlist; ranks, cell_names, and
% cell_inds are all preset to their respective values based on the Tasic
% et al, 2018 dataset used by Mezias et al, 2020, but for novel datasets
% these should be supplied by the user together. cell_names is a 1 x
% n_types cell array of descriptive names of the cell types, cell_inds is a
% 1 x n_types numeric array of the indices of the columns in
% outstruct(ind).corrB corresponding to the types in cell_names, and ranks
% is a 1 x n_types numeric array of the expected layer ordering of these
% types (e.g. a value of "1" indicates that this cell type is expected to
% be closest to the cortical surface, with ties allowed).

% Set default values if not supplied by user
if nargin < 6
    directory = [cd filesep 'MatFiles'];
    if nargin < 3
        ranks = [1 2 3 3 4 4 4];
        cell_inds = 9:15;
        cell_names = {'L2n3','L4','L5IT','L5PT','L6CT','L6IT','L6b'};
    elseif nargin < 5
        error(sprintf('Error. \nUser must supply ranks, cell names, and cell indices together or not at all.'));
    end
end

% Load dependencies
load([directory filesep 'tau_calc_dependencies.mat'],'GENGDmod','nonzerovox',...
    'structList', 'structIndex','neomask3d','neobounds3d');

% Initialize
metric_struct = struct('paramlist',{idx},'tau',...
    {zeros(length(idx),1)});
varlist = who;
varlist = [varlist;{'varlist'};{'m'}];


% Load cell count values per voxel, insert at the appropriate CCF atlas
% coordinates, store in a struct (pervoxglut)
predvecs = outstruct(idx).corrB;
newVoxMap = zeros(size(GENGDmod));
pervoxglut = struct;
for i = 1:length(cell_inds)
    pervoxglut.(cell_names{i}) = newVoxMap;
end
glutnames = cell_names;
for i = 1:length(structList)
    [~,voxinds] = ismember(structIndex{i},nonzerovox);
    allvox = nonzerovox(voxinds);
    curvox = predvecs(voxinds,9:15);
    for k = 1:length(glutnames)
        curglut = squeeze(curvox(:,k));
        pervoxglut.(glutnames{k})(allvox) = curglut;
        clear curglut
    end
    clear allvox curvox voxinds
end

% Calculate Kendall's tau for the whole brain and per slice
[d1,d2,d3] = find(neomask3d);
d1inds = unique(d1);
maplocs = d1inds;
for k = 1:length(glutnames)
    curmap = pervoxglut.(glutnames{k});
    for j = 10:length(maplocs)
        curloc = maplocs(j);
        ctxsurf = squeeze(neobounds3d(curloc,:,:));
        curneomask = squeeze(neomask3d(curloc,:,:));
        [ctxsurfx,ctxsurfy] = find(ctxsurf);
        curslc = squeeze(curmap(curloc,:,:));
        mnzs = min(nonzeros(curslc(:)));
        if isempty(mnzs)
            curslc(curslc==0) = 0;
        else
            curslc(curslc==0) = mnzs;
        end
        grayslc = mat2gray(curslc);
        binslc_whole = logical(grayslc);
        paleoinds = find(~curneomask);
        neoinds = find(curneomask);
        binslc = binslc_whole;
        binslc(paleoinds) = 0;
        skeleton = bwskel(binslc);
        [skelx,skely] = find(skeleton);
        for s = 1:length(find(skeleton))
            curx = skelx(s);
            cury = skely(s);
            dists(s) = min(sqrt((ctxsurfx-curx).^2 + (ctxsurfy-cury).^2));
        end
        if isempty(skelx) || isempty(curx) || isempty(dists)
            mndists(j-9,k) = -1;
        elseif mean(curslc(paleoinds))>mean(curslc(neoinds))
            mndists(j-9,k) = -1;
        else
            mndists(j-9,k) = mean(dists(:));
        end
        clear dists cury curx
    end
end

lranks = 0;
trueranks = 0;
for j = 1:size(mndists,1)
    curmns = mndists(j,:);
    nosig = (curmns==-1);
    n_nosig(j) = length(find(nosig));
    curmns(nosig) = [];
    [~,rankinds] = sort(curmns,'ascend');
    curranks = ranks;
    curranks(nosig) = [];
    sliceranks = curranks;
    trueranks = [trueranks curranks];
    curranks = curranks(rankinds);
    lranks = [lranks curranks];
    if isempty(curranks)
        tau_per = 0;
    else
        [tau_per,~] = corr(curranks.', sliceranks.', 'type', 'Kendall');
    end
    tau_slice(j) = (1-(n_nosig(j)/7))*tau_per;
    sqerrvec(j) = (1 - tau_per);
end

metric_struct.tau_perslice = mean(tau_slice);
lranks = lranks.';
lranks = lranks(:);
trueranks = trueranks.';
trueranks = trueranks(:);
if sum(lranks) == 0
    tau = 0;
    pval = 1;
else
    [tau,pval] = corr(lranks,trueranks,'type','Kendall');
    tau = (1-(sum(n_nosig)/(length(ranks)*length(10:length(maplocs)))))*tau;
end
metric_struct.tau = tau;
metric_struct.pval = pval;
rss = sum(sqerrvec);
N = length(maplocs);
K = outstruct.nGen;
K = K / 3855;
K = K * 100;
metric_struct.bic_pct = N*log(rss/N)+log(N)*(K);

end


