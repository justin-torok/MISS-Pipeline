function Figure_5ab_exinhplots(outstruct,idx,savenclose)

if nargin < 3
    savenclose = 0;
end

reglabs = {'Amg','Cer','Sub','Hip','Hyp','Neo','Med','Mid','Olf','Pal','Pns','Str','Tha';...
            1:11,12:22,23:25,26:36,37:56,57:94,95:119,120:140,141:148,149:156,157:169,170:177,178:212};

glutsum = sum(outstruct(idx).Bsums(1:213,[8:16 20:21]),2) + sum(outstruct(idx).Bsums(214:end,[8:16 20:21]),2);
gabasum = sum(outstruct(idx).Bsums(1:213,[1:7 17:19]),2) + sum(outstruct(idx).Bsums(214:end,[1:7 17:19]),2);
neuronratios = glutsum ./ (glutsum + gabasum); neuronratios(12) = [];

rng('default');
xpos = zeros(212,1);
g = xpos;
xposscatter = @(x,y) 0.2 * (2*rand(length(x),1) - 1) + y;
for j = 1:size(reglabs,2)
    curinds = reglabs{2,j};
    g(curinds) = j;
    xpos(curinds) = xposscatter(curinds,j);
end

cmap = hsv(size(reglabs,2));
f3 = figure('Position',[0 0 400 1000]); hold on;
v = gscatter(xpos,neuronratios,g,cmap,[],20,'off');
b = boxplot(neuronratios,g,'Colors',cmap,'Symbol','');
set(b,{'linew'},{3})
plot([0,14],[0.5,0.5],'k--','LineWidth',2);
set(gca, 'XTick', 1:13, 'XTickLabel', reglabs(1,:));
set(gca, 'YTick', [0,0.25,0.5,0.75,1])
ylabel('');
xlabel('');
set(gca,'TickLength',[0 0])
title({'Glutamatergic Fraction';'Across Regions'});
set(gca, 'FontSize', 24);
view([90 -90])
if savenclose
    print('exinhneuronratio_barplots','-dtiffn')
    close
end

cmap = lines(7);

figure('Position',[0 0 800 400]); hold on;
hg = histfit(neuronratios,10,'kernel'); hold on;
hg(1).FaceColor = cmap(4,:);
hg(1).EdgeColor = 'w';
hg(1).FaceAlpha = 0.5;
hg(2).LineWidth = 2.5;
set(gca,'FontSize',40);
set(gca,'XTick',[0.2 0.5 0.8]);
set(gca,'YTick',[0 20 40]);
title('Regional Glutamatergic Fractions','FontSize',40);
xlabel('Proportion Glutamatergic Neurons','FontSize',40);
box on;
ylim([0 50]);
xlim([0 1]);
if savenclose
    print('exinhneuronratio_hist','-dtiffn')
    close
end


end