function Figure_6a_MISS_correlations_comp(outstruct,idx,savenclose,directory)

if nargin < 4
    directory = [cd filesep 'MatFiles'];
    if nargin < 3
        savenclose = 0;
    end
end

load([directory filesep 'Zeisel_coronal_geneset.mat'],'Zeisel_gene_names','repcells','repcellinds');
load([directory filesep 'Zeisel_Inputs.mat'],'genevct','voxvgene','gene_names');

unitygenes = ismember(gene_names,Zeisel_gene_names);
E_corr = voxvgene(:,unitygenes);
C_corr = genevct(unitygenes,:);
E_corr = E_corr.';
nvox = size(E_corr,2);
logE = log2(E_corr + 1); % Zeisel et al initial log2 normalization
mn_logE = mean(logE,2); std_logE = std(logE,[],2);
Z_logE = (logE - repmat(mn_logE,1,nvox)) ./ repmat(std_logE,1,nvox);
Enans = isnan(Z_logE);
ntypes = size(C_corr,2);
mn_C = mean(C_corr,2); 
Z_C = (C_corr - repmat(mn_C,1,ntypes)) ./ repmat(mn_C,1,ntypes);
Cnans = isnan(Z_C);
allnans = logical(Cnans(:,1) + Enans(:,1)); 
Bcorr = corr(Z_C(~allnans,:),Z_logE(~allnans,:));
Bcorr = Bcorr.';

Bdense = outstruct(idx).corrB;

% newcnames = classkey(repcellinds);
% check = strcmp(repcells,newcnames);

Bdense = Bdense(:,repcellinds);
Bcorr = Bcorr(:,repcellinds);
% Bcorr(Bcorr<0) = 0;

for i = 1:size(Bdense,2)
    curdense = Bdense(:,i);
    curcorr = Bcorr(:,i);
    densetest = curdense > 0;
    corrtest = curcorr > 0;
    sumtest = densetest + corrtest;
    alltest = sumtest > 1;
    curdense = curdense(alltest);
    curcorr = curcorr(alltest);
    rval(i) = corr(curdense,curcorr,'Type','Pearson');
end

map_rvals = rval;

ctypes.Zeisel = {[143:149 150:165 35:36 21:25 45:54 81:92 9 13 18];...%Glut
    [166:186 26:33 93:106 55:63 122:126 112:117 130:131 10:11 14];...%GABA
    [16 17 118 119 107 110 111];...%Oligo
    [138 139 7 8 6 4 5 3 2];...%Astro
    [137 108 109];...%Micro/Macro
    [19:20 41:44 64:69 70:80];...%Neuromod & Neuropep
    [1 12 15 34 37:40 120:121 127:129 132:136 140:142 187:191]}; %Other

resort_ctypes = [1 2 6 3 4 5 7];
% ctypes.Zeisel = [ctypes.Zeisel;{[19:20 41:44 64:69 70:80];...
%     [1 12 15 34 37:40 120:121 127:129 132:136 140:142 187:191]}];
% ctypes.Zeisel{1}(8) = [];
% ctypes.Zeisel{3}([6 7]) = [];
% ctypes.Zeisel{5}([1 5]) = [];
ctypes.Zeisel = ctypes.Zeisel(resort_ctypes);
cmap = lines(7);
cmap = cmap(resort_ctypes,:);
cnames = {'Glut.','GABA','Nmd. & Pep.','Oligo.','Astro.','Micro.','Other'};

rng('default');
xpos = zeros(length(repcells),1);
g = xpos;
xposscatter = @(x,y) 0.2 * (2*rand(length(x),1) - 1) + y;
for j = 1:length(ctypes.Zeisel)
    curinds = ctypes.Zeisel{j};
    g(curinds) = j;
    xpos(curinds) = xposscatter(curinds,j);
end

f3 = figure('Position',[0 0 1500 400]); hold on;
v = gscatter(xpos,map_rvals,g,cmap,[],20,'off');
b = boxplot(map_rvals,g,'Colors',cmap,'Symbol','');
set(b,{'linew'},{3})
plot([0 length(cnames)+1],[median(map_rvals) median(map_rvals)],'k--','LineWidth',2);
set(gca,'XTick',1:13,'XTickLabel',cnames);
set(gca,'YTick',[0,0.25,0.5,0.75,1])
ylabel('');
xlabel('');
set(gca,'TickLength',[0 0])
title('Correlation vs. MISS Maps','FontSize',24);
set(gca, 'FontSize', 24);
% view([90 -90])
if savenclose
    print('Figure_6a_barplots','-dtiffn')
    close
end


% avgrval = mean(rval);
% medrval = median(rval);
% % figure; plot(rval);
% figure; histfit(rval,12,'kernel');
% yticks([0 10 20 30 40]);
% xticks([-0.2 0 0.2 0.4 0.6 0.8]);
% set(gca,'FontSize',18);
% xlim([-0.25 1]);
% text(-0.2,25,{sprintf('Med. R = %.2f',medrval);sprintf('Avg. R = %.2f',avgrval)},'FontSize',18);
% ylabel('Frequency','FontSize',22);
% xlabel('R-Value','FontSize',22);
% title('Correlation Vs. MISS Maps','FontSize',22);
end
