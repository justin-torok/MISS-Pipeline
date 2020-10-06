function [mistruct] = Mutual_Information_Calc(outstruct,idx,onttype,directory)

if nargin < 4
    directory = [cd filesep 'MatFiles'];
    if nargin < 3
        onttype = 'fore';
    end
end
ontstr = ['ontologystruct_' onttype];
ontstruct = load([directory filesep ontstr '.mat'], ontstr);
if ~isempty(outstruct)
    [~, cell_types] = Voxel_To_Region_Bilateral(outstruct(idx).corrB,directory);
else
    ISHdat = load('PresetInputs.mat','regvgene');
    ISHdat = ISHdat.regvgene;
    [~, cell_types] = Voxel_To_Region_Bilateral(ISHdat,directory);
end
mistruct = struct;
mistruct.onttype = onttype;

% rng('default');
iters = 10000;

nclust = ontstruct.(ontstr).ABI.nclust(1:end);
clustmat = ontstruct.(ontstr).ABI.clusters;
miindex = zeros(length(nclust),iters);
miT_adj = zeros(1,length(nclust));
stds = zeros(1,length(nclust));
% metric1 = zeros(1,length(nclust));

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
    miT = mutInfo(Tclust,vec);
    miTs(n) = miT;
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
        miindex(n,k) = mutInfo(randvec,vec);
        clear randvec randszvec
    end
    stds(n) = std(miindex(n,:));
    miT_adj(n) = (miT - mean(miindex(n,:)))/(mean([entropy(vec),entropy(Tclust)]) - mean(miindex(n,:)));
%     metric1(n) = (miT - mean(miindex(n,:))) / stds(n);
end
mistruct.nclust = nclust;
mistruct.MutualInfoT = miTs;
mistruct.MutualInfoRandom = miindex;
mistruct.AdjustedMutualInfo = miT_adj;
% mistruct.StdAboveMean = metric1;
mistruct.Tmat = Tmat;

end


