function [randstruct] = Rand_Index_Calc(outstruct,idx,onttype,directory)

if nargin < 4
    directory = [cd filesep 'MatFiles'];
end
ontstr = ['ontologystruct_' onttype];
ontstruct = load([directory filesep ontstr '.mat'], ontstr);
[~, cell_types] = Voxel_To_Region_Bilateral(outstruct(idx).corrB);
randstruct = struct;
randstruct.onttype = onttype;

% rng('default');
iters = 10000;

nclust = ontstruct.(ontstr).ABI.nclust(1:end);
clustmat = ontstruct.(ontstr).ABI.clusters;
randindex = zeros(length(nclust),iters);
randT_adj = zeros(1,length(nclust));
stds = zeros(1,length(nclust));
metric1 = zeros(1,length(nclust));

amg = 1:11; sub = 23:25; hip = 26:36; hyp = 37:56; neo = 57:94; 
olf = 141:148; pal = 149:156; str = 170:177; tha = 178:212;
cer = 12:22; med = 95:119; mid = 120:140; pns = 157:169;

X = cell_types;
X(12,:) = [];

if strcmp(onttype,'fore')
    X([cer med mid pns],:) = [];
elseif strcmp(onttype,'hind')
    X([amg sub hip hyp neo olf pal str tha],:) = [];
end
Y = pdist(X);
Z = linkage(Y, 'ward');

for n = 1:length(nclust)
    Tclust = cluster(Z, 'maxclust', nclust(n));
    Tmat(:,n) = Tclust;
    vec = clustmat(:,n); % ground truth ontological clustering
    ubmat = zeros(length(vec));
    uTmat = zeros(length(Tclust));
    for i = 1:length(vec)
        for j = i:length(vec)
            if j > i
                ubmat(i,j) = 2*(vec(i)==vec(j))-1;
                uTmat(i,j) = 2*(Tclust(i)==Tclust(j))-1;
            end
        end
    end
    triu_indmat = ones(length(Tclust)) - tril(ones(length(Tclust)),0);
    triu_ind = logical(triu_indmat(:));
    ubvec = ubmat(:); ubvec = ubvec(triu_ind);
    uTvec = uTmat(:); uTvec = uTvec(triu_ind);
    uvec_true = ubvec .* uTvec;
    randT = sum(uvec_true==1)/length(uvec_true);
    randTs(n) = randT;
    trueszvec = zeros(1,nclust(n));
    for n2 = 1:nclust(n)
        trueszvec(n2) = sum(vec == n2);
    end
    trueszvec = sort(trueszvec); % make sure the largest group is last
    for k = 1:iters
        randszvec = zeros(1,nclust(n));
        randvec = [];
        for n2 = 1:(nclust(n))
            if n2 == nclust(n)
                randszvec(n2) = length(Tclust) - sum(randszvec);
            else
                randszvec(n2) = randi([round(0.8*trueszvec(n2)),round(1.2*trueszvec(n2))]);
            end
            randvec = [randvec, n2*ones(1,randszvec(n2))];
        end
        randinds = randperm(length(randvec));
        randvec = randvec(randinds);
        uamat = zeros(length(randvec));
        for i = 1:length(randvec)
            for j = i:length(randvec)
                if j > i
                    uamat(i,j) = 2*(randvec(i)==randvec(j))-1;
                end
            end
        end
        uavec = uamat(:); uavec = uavec(triu_ind);
        uvec = uavec .* ubvec;
        randindex(n,k) = sum(uvec==1)/length(uvec);
        clear randvec randszvec
    end
    stds(n) = std(randindex(n,:));
    randT_adj(n) = (randT - mean(randindex(n,:)))/(1 - mean(randindex(n,:)));
    metric1(n) = (randT - mean(randindex(n,:))) / stds(n);
end
randstruct.nclust = nclust;
randstruct.RandIndexT = randTs;
randstruct.RandIndexRandom = randindex;
randstruct.AdjustedRand = randT_adj;
randstruct.StdAboveMean = metric1;
randstruct.Tmat = Tmat;

end


