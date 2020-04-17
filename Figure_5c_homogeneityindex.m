function Figure_5c_homogeneityindex(outstruct,idx,savenclose)
if nargin < 3
    savenclose = 0;
end
    
regtotals = outstruct(idx).Bsums;
normtotals = zeros(426,25);
regentropy = zeros(1,426);
entropycalc = @(x) -dot(log(x+eps),x);
for i = 1:426
   normtotals(i,:) = regtotals(i,:)/sum(regtotals(i,:)); 
   regentropy(i) = entropycalc(normtotals(i,:));
end
regentropy(isnan(regentropy)) = 0;
datinput = ones(426,1)./regentropy.';
datinput([12 12+213]) = [];

newinput = (datinput(1:212) + datinput(213:end)) / 2;
reglabs = {'Amg','Cer','Sub','Hip','Hyp','Neo','Med','Mid','Olf','Pal','Pns','Str','Tha';...
            1:11,12:22,23:25,26:36,37:56,57:94,95:119,120:140,141:148,149:156,157:169,170:177,178:212};
plotlabs = reglabs(1,:);
rng('default');
xpos = zeros(212,1);
g = xpos;
xposscatter = @(x,y) 0.2 * (2*rand(length(x),1) - 1) + y;
for j = 1:length(plotlabs)
    curinds = reglabs{2,j};
    g(curinds) = j;
    xpos(curinds) = xposscatter(curinds,j);
end

cmap = hsv(13);
f3 = figure('Position',[0 0 500 1000]); hold on; box on;
v = gscatter(newinput,xpos,g,cmap,[],20,'off');
b = boxplot(newinput,g,'orientation','horizontal','Colors',cmap,'Symbol','');
set(b,{'linew'},{3})
plot([median(newinput),median(newinput)],[0,14],'k--','LineWidth',2);
set(gca, 'YTick', 1:13, 'YTickLabel', reglabs(1,:));
set(gca, 'XTick', [0.5 2 3.5])
ylabel('');
xlabel('');
xlim([0 max(datinput)]);
set(gca,'TickLength',[0 0])
title({sprintf('Homogeneity Index \n Across Regions')});
set(gca, 'FontSize', 30,'LineWidth',0.75);

if savenclose
    print('Figure_5c_homogeneityindex','-dtiffn');
    close
end
end