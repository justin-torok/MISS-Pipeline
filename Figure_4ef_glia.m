function Figure_4ef_glia(outstruct,idx,savenclose,directory)
% Designed to handle an outstruct of dimension length(ngenelist) that has a
% single field "bvals." idx is the index of interest that is within the
% range of ngenelist
if nargin < 4
    directory = [cd filesep 'MatFiles'];
    if nargin < 3
        savenclose = 0;
    end
end

load([directory filesep 'Tasic_Inputs.mat'],'listBmap');

neoinds = 57:94; 
foreinds = [1:11,23:56,141:156,170:212];
wbinds = [12:22 95:140 157:169];
indcell = {wbinds,foreinds,neoinds};
indcolor = {'b','m','g'};
indmark = {'*','o','d','s','^','p'};

LinRcalc = @(x,y) 2*corr(x,y)*std(x)*std(y)/(std(x)^2 + std(y)^2 + (mean(x) - mean(y))^2);

bigmicro = 0; datmicro = 0;
listBmap = listBmap(:);
for i = 1:212
    voxels(i) = sum(listBmap==i);
end
reglabs = {'Amg','Cer','Sub','Hip','Hyp','Neo','Med','Mid','Olf','Pal','Pns','Str','Tha';...
            1:11,12:22,23:25,26:36,37:56,57:94,95:119,120:140,141:148,149:156,157:169,170:177,178:212};
astrosum = sum(outstruct(idx).Bsums(1:213,22),2) + sum(outstruct(idx).Bsums(214:end,22),2);
microsum = sum(outstruct(idx).Bsums(1:213,23),2) + sum(outstruct(idx).Bsums(214:end,23),2);
oligosum = sum(outstruct(idx).Bsums(1:213,24),2) + sum(outstruct(idx).Bsums(214:end,24),2);
for i = 1:size(reglabs,2)
    curinds = reglabs{2,i};
    astroplot(i) = 125*sum(astrosum(curinds))/sum(voxels(curinds));
    microplot(i) = 125*sum(microsum(curinds))/sum(voxels(curinds));
    oligoplot(i) = 125*sum(oligosum(curinds))/sum(voxels(curinds));
end
gliaplot = [astroplot; microplot; oligoplot]; gliaplot = gliaplot.';


f1 = figure('Units','inch','Position',[0 0 16 4]); hold on;

for i = 1:length(indcell)
    load([directory filesep 'keller_micro_totals_reorder.mat'],'keller_micro_totals_reorder');
    load([directory filesep 'keller_micro_listB.mat'],'keller_micro_listB');
    testinds = indcell{i};    
    sumcts_wb = outstruct(idx).Bsums(1:213,:) + outstruct(idx).Bsums(214:end,:);
    sumcts_wb(12,:) = [];
    sumcts_micro_wb = sumcts_wb(:,23);

    hh_micro = [];
    infer_micro = [];
    for j = 1:size(keller_micro_listB)
        hh = keller_micro_listB(j,:);
        hhinds = hh(hh > 0);
        if all(ismember(hhinds,testinds))
            hh_micro = [hh_micro, sum(keller_micro_totals_reorder(hhinds))];
            infer_micro = [infer_micro, sum(sumcts_micro_wb(hhinds))];
        end
    end
    
    set(gca,'FontSize',18);
    subplot(1,3,1); hold on;
    scatter(infer_micro,hh_micro,30,[indcolor{i} indmark{5}],'filled'); hold on;
    ylabel('Empirical Counts','FontSize',20)
    xlabel('Inferred Counts','FontSize',20)
    
    bigmicro = [bigmicro; infer_micro.'];
    datmicro = [datmicro; hh_micro.'];

end

bigmicro(1) = []; datmicro(1) = [];

subplot(1,3,1); 
hold on; box on;
title('Microglia','FontSize',20)
plot([1 max(cat(1,bigmicro,datmicro))],[1 max(cat(1,bigmicro,datmicro))],'k-'); hold on;
txt = {sprintf('R_C = %.2f',LinRcalc(bigmicro,datmicro));...
    sprintf('R = %.2f',corr(bigmicro,datmicro))};
text(2.2*10^5,0.5*10^5,txt,'FontSize',15);
maxlab = max(cat(1,bigmicro,datmicro));
set(gca,'YTick',[0 0.5*maxlab maxlab]);
set(gca,'YTickLabels',{'0',sprintf('%.1f',0.5*maxlab/10^5),'x10^5'});
set(gca,'XTick',[0 0.5*maxlab maxlab]);
set(gca,'XTickLabels',{'0',sprintf('%.1f',0.5*maxlab/10^5),'x10^5'});
p = polyfit(bigmicro,datmicro,1);
x_maxy = (maxlab - p(2))/p(1);
if x_maxy < maxlab
    plot([0,x_maxy],[p(2),(p(2) + x_maxy*p(1))],'Color',[0.9 0.3 0.1]); hold on;
else
    plot([0,maxlab],[p(2),(p(2) + maxlab*p(1))],'Color',[0.9 0.3 0.1]); hold on;
end
xlim([0 max(cat(1,bigmicro,datmicro))]);
ylim([0 max(cat(1,bigmicro,datmicro))]);

subplot(1,3,[2,3])
bar(gliaplot);
set(gca,'XTickLabel',reglabs(1,:));
set(gca,'YTick',[0,10000,20000]);
ylabel('Cell Density (cells/mm^3)');
set(gca,'TickLength',[0 0]);
legend('Astro','Micro','Oligo','Location','northeast');
set(gca,'FontSize',18);
title('Glial Density Per Major Region Group','FontSize',20);

if savenclose
    print('Figure_4ef_gliaplots','-dtiffn');
    close
end
end