function T_z = Figure_7a_dendrograms_Zeisel(outstruct,idx,onttype,savenclose,directory)

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
T = cluster(Z,'MaxClust',5);
T_z = T;
optimleaford = optimalleaforder(Z,Y);
figure('Units', 'inch', 'Position', [0 0 15 6]);
[D1,~,outperm] = dendrogram(Z,0,'Orientation','top','reorder',optimleaford); hold on;
set(D1,'LineWidth',2);
set(D1,'Color','k');
set(gca,'XTick',[]);
set(gca,'YTick',[]);
f1 = find(T==2); f2 = find(T==1); f3 = find(T==3); f4 = find(T==4); f5 = find(T==5);
is1 = ismember(outperm,f1); f = find(is1); xlimits(1) = f(end); 
is2 = ismember(outperm,f2); f = find(is2); xlimits(2) = f(end); 
is3 = ismember(outperm,f3); f = find(is3); xlimits(3) = f(end); 
is4 = ismember(outperm,f4); f = find(is4); xlimits(4) = f(end);
is5 = ismember(outperm,f5); f = find(is5); xlimits(5) = f(end);
[xlimits,xlimidx] = sort(xlimits);
xlimdiffs = diff(xlimits);
xlimdiffs = [xlimits(1) xlimdiffs];
% splits = [is1.' is2.' is3.' is4.' is5.'];
% splits = splits * 3000;
areacols = lines(length(unique(T)));
R1 = rectangle('Position',[0 0 xlimits(1)+0.5 4000],'FaceColor',[areacols(1,:) 0.6],'EdgeColor','w'); hold on;
R2 = rectangle('Position',[xlimits(1)+0.5 0 xlimdiffs(2) 4000],'FaceColor',[areacols(5,:) 0.6],'EdgeColor','w'); hold on;
R3 = rectangle('Position',[xlimits(2)+0.5 0 xlimdiffs(3) 4000],'FaceColor',[areacols(3,:) 0.6],'EdgeColor','w'); hold on;
R4 = rectangle('Position',[xlimits(3)+0.5 0 xlimdiffs(4) 4000],'FaceColor',[areacols(4,:) 0.6],'EdgeColor','w'); hold on;
R5 = rectangle('Position',[xlimits(4)+0.5 0 xlimdiffs(5) 4000],'FaceColor',[areacols(2,:) 0.6],'EdgeColor','w'); hold on;
% A1 = area([(0:size(X,1)+1).' (0:size(X,1)+1).' (0:size(X,1)+1).' (0:size(X,1)+1).' (0:size(X,1)+1).'],...
% [0 0 0 3000 0;splits;0 0 0 0 3000],'FaceAlpha',0.6,'EdgeAlpha',0.1,'LineWidth',2.5); hold on;
% A1 = area([0 0 0 3000 0;splits;0 0 0 0 3000],'FaceAlpha',0.6,'EdgeAlpha',0.1,'LineWidth',2.5); hold on;
% areacols = lines(5);
% set(A1,{'facecolor'},{areacols(1,:);areacols(2,:);areacols(3,:);areacols(4,:);areacols(5,:)});
% set(A1,'EdgeColor','w');
ax = gca;
set(ax,'YColor',[1 1 1]);
xlim([0 size(X,1)+1]);
ylim([0 4000]);
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
T = cluster(Z,'MaxClust',5);
optimleaford = optimalleaforder(Z,Y);
figure('Units', 'inch', 'Position', [0 0 15 6]);
[D2,~,outperm] = dendrogram(Z,size(clusts,1),'reorder',optimleaford); hold on;
set(D2,'LineWidth',2);
set(D2,'Color','k');
set(gca,'XTick',[]);
set(gca,'YTick',[]);
f1 = find(T==2); f2 = find(T==1); f3 = find(T==3); f4 = find(T==4); f5 = find(T==5);
is1 = ismember(outperm,f1); f = find(is1); xlimits(1) = f(end); 
is2 = ismember(outperm,f2); f = find(is2); xlimits(2) = f(end); 
is3 = ismember(outperm,f3); f = find(is3); xlimits(3) = f(end); 
is4 = ismember(outperm,f4); f = find(is4); xlimits(4) = f(end);
is5 = ismember(outperm,f5); f = find(is5); xlimits(5) = f(end);
[xlimits,xlimidx] = sort(xlimits);
xlimdiffs = diff(xlimits);
xlimdiffs = [xlimits(1) xlimdiffs];
areacols = lines(length(unique(T)));
R1 = rectangle('Position',[0 0 xlimits(1)+0.5 4000],'FaceColor',[areacols(2,:) 0.6],'EdgeColor','w'); hold on;
R2 = rectangle('Position',[xlimits(1)+0.5 0 xlimdiffs(2) 4000],'FaceColor',[areacols(4,:) 0.6],'EdgeColor','w'); hold on;
R3 = rectangle('Position',[xlimits(2)+0.5 0 xlimdiffs(3) 4000],'FaceColor',[areacols(1,:) 0.6],'EdgeColor','w'); hold on;
R4 = rectangle('Position',[xlimits(3)+0.5 0 xlimdiffs(4) 4000],'FaceColor',[areacols(3,:) 0.6],'EdgeColor','w'); hold on;
R5 = rectangle('Position',[xlimits(4)+0.5 0 xlimdiffs(5) 4000],'FaceColor',[areacols(5,:) 0.6],'EdgeColor','w'); hold on;
% splits = [is1.' is2.' is3.' is4.' is5.'];
% splits = splits * 5.5;
% A2 = area([(0:size(clusts,1)+1).' (0:size(clusts,1)+1).'],[5.5 0;splits;0 5.5],'FaceAlpha',0.1,'EdgeAlpha',0.1,'LineWidth',2.5); hold on;
% area_colors = lines(5);
% area_colors = mat2cell(area_colors);
% set(A2,{'facecolor'},area_colors);
% set(A2,'EdgeColor','w');
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