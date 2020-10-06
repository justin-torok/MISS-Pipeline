function Figure_6a_dendrograms(outstruct,idx,onttype,savenclose,directory)

if nargin < 5
    directory = [cd filesep 'MatFiles'];
    if nargin < 4
        savenclose = 0;
        if nargin < 3
            onttype = 'fore';
        end
    end
end

ontstr = ['ontologystruct_' onttype];
ontstruct = load([directory filesep ontstr '.mat'], ontstr);
clusts = ontstruct.(ontstr).ABI.clusters;

[~, cell_types] = Voxel_To_Region_Bilateral(outstruct(idx).corrB, directory);
X = cell_types; X(12,:) = [];

M(1:11) = 1; M(12:23) = 2; M(24:26) = 3; M(27:37) = 4; M(38:57) = 5;
M(58:95) = 6; M(96:120) = 7; M(121:141) = 8; M(142:149) = 9;
M(150:157) = 10; M(158:170) = 11; M(171:178) = 12; M(179:213) = 13;
M(12) = [];
cmap = colormap(hsv(length(unique(M))));

amg = 1:11; sub = 23:25; hip = 26:36; hyp = 37:56; neo = 57:94; 
olf = 141:148; pal = 149:156; str = 170:177; tha = 178:212;
cer = 12:22; med = 95:119; mid = 120:140; pns = 157:169;

if strcmp(onttype,'fore')
    M([cer med mid pns]) = [];
    uniM = unique(M);
    for j = 1:length(unique(M))
        M(M==uniM(j)) = j;
    end
    X([cer med mid pns],:) = [];
    cmap = cmap([1 3 4 5 6 9 10 12 13],:);
elseif strcmp(onttype,'hind')
    M([amg sub hip hyp neo olf pal str tha]) = [];
    uniM = unique(M);
    for j = 1:length(unique(M))
        M(M==uniM(j)) = j;
    end
    X([amg sub hip hyp neo olf pal str tha],:) = [];
    cmap = cmap([2 7 8 11],:);
end

Y = pdist(X);
Z = linkage(Y, 'ward');
T = cluster(Z,'MaxClust',2);
optimleaford = optimalleaforder(Z,Y);
figure('Units', 'inch', 'Position', [0 0 15 6]);
[D1,~,outperm] = dendrogram(Z,0,'Orientation','top','reorder',optimleaford); hold on;
set(D1,'LineWidth',2);
set(D1,'Color','k');
set(gca,'XTick',[]);
set(gca,'YTick',[]);
f1 = find(T==2); f2 = find(T==1);
is1 = ismember(outperm,f1);
is2 = ismember(outperm,f2);
splits = [is1.' is2.'];
splits = splits * 3000;
A1 = area([(0:size(X,1)+1).' (0:size(X,1)+1).'],[3000 0;splits;0 3000],'FaceAlpha',0.1,'EdgeAlpha',0.1,'LineWidth',2.5); hold on;
set(A1,{'facecolor'},{'g';'m'});
set(A1,'EdgeColor','w');
ax = gca;
set(ax,'YColor',[1 1 1]);
xlim([0 size(X,1)+1]);
ylim([0 3000]);
title('Cell-Type Based Clustering','FontSize',30);
if savenclose
    print('dendogramcelltype','-dtiffn');
    close
end
outperm = outperm.';
MM = M(outperm);
figure('Units','inch','Position',[0 0 14.9 0.25]);
imagesc(MM); colormap(cmap);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
if savenclose
    print('dendogramcelltype_colorbar','-dtiffn');
    close
end

rng('default');
pdistmat = zeros(size(clusts,1));
for i = 1:size(pdistmat,1)
    for j = 1:size(pdistmat,2)
        if j > i
            pdistmat(i,j) = 0.8*rand;
            for k = 1:size(clusts,2)
                clust = clusts(:,k);
                if clust(i) ~= clust(j)
                    pdistmat(i,j) = pdistmat(i,j) + 1;
                end
            end
        end
    end
end

pdistmat = pdistmat + pdistmat.';
Y = squareform(pdistmat);
Z = linkage(Y,'average');
T = cluster(Z,'MaxClust',2);
optimleaford = optimalleaforder(Z,Y);
figure('Units', 'inch', 'Position', [0 0 15 6]);
[D2,~,outperm] = dendrogram(Z,size(clusts,1),'reorder',optimleaford); hold on;
set(D2,'LineWidth',2);
set(D2,'Color','k');
set(gca,'XTick',[]);
set(gca,'YTick',[]);
f1 = find(T==2); f2 = find(T==1);
is1 = ismember(outperm,f1);
is2 = ismember(outperm,f2);
splits = [is1.' is2.'];
splits = splits * 5.5;
A2 = area([(0:size(clusts,1)+1).' (0:size(clusts,1)+1).'],[5.5 0;splits;0 5.5],'FaceAlpha',0.1,'EdgeAlpha',0.1,'LineWidth',2.5); hold on;
set(A2,{'facecolor'},{'g';'m'});
set(A2,'EdgeColor','w');
ax = gca;
set(ax,'YColor',[1 1 1]);
xlim([0 size(clusts,1)+1]);
ylim([0 5.5]);
title('Developmental Ontology Divisions','FontSize',30);
if savenclose
    print('dendogramontology','-dtiffn');
    close
end
outperm = outperm.';
MM = M(outperm);
figure('Units','inch','Position',[0 0 14.9 0.25]);
imagesc(MM); colormap(cmap);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
if savenclose
    print('dendogramontology_colorbar','-dtiffn');
    close
end

end