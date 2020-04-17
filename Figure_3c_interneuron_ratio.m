function Figure_3c_interneuron_ratio(outstruct,idx,savenclose,directory)

if nargin < 4
    directory = [cd filesep 'MatFiles'];
    if nargin < 3
        savenclose = 0;
    end
end

load([directory filesep 'PresetInputs.mat'],'listBmap');

listBmap = listBmap(:);
for i = 1:212
    voxels(i) = sum(listBmap==i);
end

% Interneuron Panel

pvalbsum = sum(outstruct(idx).Bsums(1:213,3),2) + sum(outstruct(idx).Bsums(214:end,3),2);
sstsum = sum(outstruct(idx).Bsums(1:213,6),2) + sum(outstruct(idx).Bsums(214:end,6),2);
vipsum = sum(outstruct(idx).Bsums(1:213,7),2) + sum(outstruct(idx).Bsums(214:end,7),2);
intersum = pvalbsum + sstsum + vipsum;

zerosum = (intersum==0);
dummysum = zeros(length(intersum),1);
dummysum(zerosum) = 1;

pvalbprop = pvalbsum ./ intersum;
sstprop = sstsum ./ intersum;
vipprop = vipsum ./ intersum;

interprop = [pvalbprop sstprop vipprop dummysum];
intersum(zerosum) = min(intersum)/2;
interprop(12,:) = []; intersum(12,:) = [];

reglabs = {'Amg','Cer','Sub','Hip','Hyp','Neo','Med','Mid','Olf','Pal','Pns','Str','Tha';...
            1:11,12:22,23:25,26:36,37:56,57:94,95:119,120:140,141:148,149:156,157:169,170:177,178:212};
plotlabs = reglabs(1,:);
for i = 1:size(reglabs,2)
    curinds = reglabs{2,i};
    newinds = curinds + i - 1;
    plotprop_i(newinds,:) = interprop(curinds,:);
    interplot(newinds) = intersum(curinds);
    xinds(i) = (newinds(end)-newinds(1)) / 2 + newinds(1);
end
cmap = flipud(pink(40));
legplot = 1:40;
figure('Position',[0 0 1000 50]);
imagesc(legplot); colormap(cmap);
set(gca,'XTick',[0, ]);
set(gca,'YTick',[]);

f1 = figure('Position',[0 0 1500 800]); hold on;
subplot(8,15,1,'Parent',f1);
subplot(8,15,1:7*15); bar(plotprop_i,'stacked','LineWidth',1); hold on;
title('Ratios of Pvalb+, Sst+, and Vip+ Interneurons','FontSize',24);
set(gca,'FontSize',20);
set(gca,'XTick',[]);
set(gca,'YTick',[0.2 0.4 0.6 0.8 1]);
ylabel('Normalized Ratio','FontSize',24);
set(gca,'TickLength',[0 0]);
ylim([0 1]);
xlim([0.5 newinds(end)+0.5]);
legend('Pvalb','Sst','Vip','None','Location','southeast');
subplot(8,15,7*15+1:8*15); 
imagesc(log(interplot+1),[min(log(intersum+1))-1 max(log(intersum+1))]); 
colormap(cmap); hold on;
set(gca,'XTick',xinds);
set(gca,'XTickLabel',plotlabs);
set(gca,'FontSize',18);
set(gca,'TickLength',[0 0]);
set(gca,'YTick',[]);

if savenclose
    print('Figure_3c_interneuron_ratios','-dtiffn');
    close
    print('Figure_3c_interneuron_ratios_heatmap','-diffn');
    close
end
end