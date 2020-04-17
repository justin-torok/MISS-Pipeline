function Figure_6d_ribrainframe(randstruct,savenclose,directory)
% Note - does not currently work for mid/hindbrain randstruct
if nargin < 3
    directory = [cd filesep 'MatFiles'];
    if nargin < 2
        savenclose = 1; % many brainframe images being rendered, better to close
    end
end
load([directory filesep 'input_struct_randindex.mat'],'input_struct')
ontstr = ['ontologystruct_' randstruct.onttype];
ontstruct = load([directory filesep ontstr '.mat'], ontstr);

foreinds = {1:11,12:14,15:25,26:45,46:83,84:91,92:99,100:107,108:142};
wholeinds = {1:11,12:23,24:26,27:37,38:57,58:95,96:120,121:141,142:149,150:157,158:170,171:178,179:213};
isfore = [1 0 1 1 1 1 0 0 1 1 0 1 1];

Tmat = randstruct.Tmat;

for i = 1:size(Tmat,2)
    nclust = ontstruct.(ontstr).ABI.nclust(i);
    T_ont = ontstruct.(ontstr).ABI.clusters(:,i);
    T_clust = Tmat(:,i);
    agreemat = zeros(nclust);
    nsumsvec = zeros(1,nclust);
    for j = 1:nclust
        for k = 1:nclust
            Taggj = (T_ont == j);
            Taggk = (T_clust == k);
            agreemat(j,k) = sum(Taggj .* Taggk);
        end
        nsumsvec(j) = sum(T_ont == j);
    end
    [~,sortind] = sort(nsumsvec,'descend');
    agreemat_ = agreemat(sortind,:);
    for j = 1:nclust
        [~,maxind] = max(agreemat_(j,:));
        agreemat_(:,maxind) = -1;
        T_clust(Tmat(:,i) == maxind) = sortind(j);
    end
    Tmat(:,i) = T_clust;
end

Tont = ontstruct.(ontstr).ABI.clusters;
newT = zeros(1,size(Tmat,2));
newOnt = zeros(1,size(Tont,2));
dex = 1;
for i = 1:length(wholeinds)
    curinds = wholeinds{i};
    if isfore(i)
        curfore = foreinds{dex};
        newT = [newT;Tmat(curfore,:)];
        newOnt = [newOnt;Tont(curfore,:)];
        dex = dex + 1;
    else
        newT = [newT;zeros(length(curinds),size(Tmat,2))];
        newOnt = [newOnt;zeros(length(curinds),size(Tont,2))];
    end
end

newT(1,:) = [];
newT = [newT;newT];
newOnt(1,:) = [];
newOnt = [newOnt;newOnt];
scalevec = ones(length(newT),1);
cmap = hsv(13);

cmap_cell = {cmap([6,13],:);cmap([6,13,8],:);cmap([6,13,8,1],:);...
    cmap([6,13,8,1,3],:);cmap([6,13,8,1,2,11,10,9,3,7,4],:)};
fnms = {'k2','k3','k4','k5','k11'};
input_struct.savenclose = savenclose;
input_struct.xfac = 0.01;
input_struct.pointsize = 10;
for i = 1:size(newT,2)
    curT = newT(:,i);
    curcol = cmap_cell{i};
    curfnm = fnms{i};
    input_struct.region_groups = curT;
    input_struct.cmap = curcol;
    input_struct.img_labels = curfnm;
    figure
    brainframe(input_struct);
end

ont_fnms = {'Dat2','Dat3','Dat4','Dat5','Dat11'};
for i = 1:size(newOnt,2)
    curOnt = newOnt(:,i);
    curcol = cmap_cell{i};
    curfnm_ont = ont_fnms{i};
    input_struct.region_groups = curOnt;
    input_struct.cmap = curcol;
    input_struct.img_labels = curfnm_ont;
    figure
    brainframe(input_struct);
end
end