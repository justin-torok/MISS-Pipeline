function Figure_4b_interneuron_scatter(outstruct,idx,savenclose,directory)
% Designed to handle an outstruct of dimension length(ngenelist) that has a
% single field "bvals." idx is the index of interest that is within the
% range of ngenelist

if nargin < 4
    directory = [cd filesep 'MatFiles'];
    if nargin < 3
        savenclose = 0;
    end
end

sampleinds = [70,92,186]; % Secondary motor area, primary visual area, dorsal lateral geniculate 
neoinds = [57:69, 71:91, 93, 94]; 
foreinds = [1:11,23:56,141:156,170:185, 187:212];
wbinds = [12:22 95:140 157:169];
indcell = {wbinds,foreinds,neoinds,sampleinds};
indcolor = {'b','m','g','r'};
indsize =  [30,30,30,100];
indmark = {'s','o','d','*','^','p'};

LinRcalc = @(x,y) 2*corr(x,y)*std(x)*std(y)/(std(x)^2 + std(y)^2 + (mean(x) - mean(y))^2);

bigpv = 0; datpv = 0;
bigsst = 0; datsst = 0;
bigvip = 0; datvip = 0;

f1 = figure('Units','inch','Position',[0 0 16 4]);
s1 = subplot(3,1,1,'Parent',f1); hold on; 
for i = 1:length(indcell)
    load([directory filesep 'kim_totals_reorder_m.mat'],'kim_totals_reorder');
    testinds = indcell{i};
    kim_totals_reorder = kim_totals_reorder(testinds,:);
    nonnaninds = zeros(size(kim_totals_reorder));
    for j = 1:size(kim_totals_reorder,2)
        nonnaninds(:,j) = ~isnan(kim_totals_reorder(:,j));
    end
    nonnaninds = logical(nonnaninds);
    kim_pv_wb = kim_totals_reorder(nonnaninds(:,1),1);
    kim_sst_wb = kim_totals_reorder(nonnaninds(:,2),2);
    kim_vip_wb = kim_totals_reorder(nonnaninds(:,3),3);

    sumcts_wb = outstruct(idx).Bsums(1:213,:) + outstruct(idx).Bsums(214:end,:);
    sumcts_wb(12,:) = [];
    sumcts_wb_mod = sumcts_wb(testinds,:);
    sumcts_pv_wb = sumcts_wb_mod(nonnaninds(:,1),3);
    sumcts_sst_wb = sumcts_wb_mod(nonnaninds(:,2),6);
    sumcts_vip_wb = sumcts_wb_mod(nonnaninds(:,3),7);
    
    cellsim_cell = {sumcts_pv_wb, sumcts_sst_wb, sumcts_vip_wb};          
    celldat_cell = {kim_pv_wb, kim_sst_wb, kim_vip_wb};
                
    bigpv = [bigpv; sumcts_pv_wb];
    datpv = [datpv; kim_pv_wb];
    bigsst = [bigsst; sumcts_sst_wb];
    datsst = [datsst; kim_sst_wb];
    bigvip = [bigvip; sumcts_vip_wb];
    datvip = [datvip; kim_vip_wb];

    subplot(1,3,1); hold on;
    set(gca,'FontSize',12);
    scatter(cellsim_cell{1},celldat_cell{1},indsize(i),[indcolor{i} indmark{1}],'filled'); hold on;
    
    subplot(1,3,2); hold on;
    set(gca,'FontSize',12);
    scatter(cellsim_cell{2},celldat_cell{2},indsize(i),[indcolor{i} indmark{2}],'filled'); hold on;
    
    subplot(1,3,3); hold on;
    set(gca,'FontSize',12);
    scatter(cellsim_cell{3},celldat_cell{3},indsize(i),[indcolor{i} indmark{3}],'filled'); hold on;

end

bigpv(1) = []; bigsst(1) = []; bigvip(1) = [];
datpv(1) = []; datsst(1) = []; datvip(1) = []; 

subplot(1,3,1); hold on; box on;
title('Pvalb','FontSize',20)
plot([1 max([bigpv;datpv])],[1 max(cat(1,bigpv,datpv))],'k-'); hold on;
txt = {sprintf('R_C = %.2f',LinRcalc(bigpv,datpv));...
    sprintf('R = %.2f',corr(bigpv,datpv))};
text(1.1*10^5,2.4*10^5,txt,'FontSize',15);
ylabel('Empirical Counts','FontSize',20);
xlabel('Inferred Counts','FontSize',20);
maxlab = max(cat(1,bigpv,datpv));
set(gca,'YTick',[0 0.5*maxlab maxlab]);
set(gca,'YTickLabels',{'0',sprintf('%.1f',0.5*maxlab/10^5),'x10^5'});
set(gca,'XTick',[0 0.5*maxlab maxlab]);
set(gca,'XTickLabels',{'0',sprintf('%.1f',0.5*maxlab/10^5),'x10^5'});
set(gca,'FontSize',18);
p = polyfit(bigpv,datpv,1);
x_maxy = (maxlab - p(2))/p(1);
if x_maxy < maxlab
    plot([0,x_maxy],[p(2),(p(2) + x_maxy*p(1))],'Color',[0.9 0.3 0.1]); hold on;
else
    plot([0,maxlab],[p(2),(p(2) + maxlab*p(1))],'Color',[0.9 0.3 0.1]); hold on;
end
xlim([0 max(cat(1,bigpv,datpv))]);
ylim([0 max(cat(1,bigpv,datpv))]);

subplot(1,3,2); hold on; hold on; box on;
title('Sst','FontSize',20)
plot([1 max(cat(1,bigsst,datsst))],[1 max(cat(1,bigsst,datsst))],'k-'); hold on;
txt = {sprintf('R_C = %.2f',LinRcalc(bigsst,datsst));...
    sprintf('R = %.2f',corr(bigsst,datsst))};
text(3.75*10^4,8.5*10^4,txt,'FontSize',15);
xlabel('Inferred Counts','FontSize',20);
maxlab = max(cat(1,bigsst,datsst));
set(gca,'YTick',[0 0.5*maxlab maxlab]);
set(gca,'YTickLabels',{'0',sprintf('%.1f',0.5*maxlab/10^4),'x10^4'});
set(gca,'XTick',[0 0.5*maxlab maxlab]);
set(gca,'XTickLabels',{'0',sprintf('%.1f',0.5*maxlab/10^4),'x10^4'});
set(gca,'FontSize',18);
p = polyfit(bigsst,datsst,1);
x_maxy = (maxlab - p(2))/p(1);
if x_maxy < maxlab
    plot([0,x_maxy],[p(2),(p(2) + x_maxy*p(1))],'Color',[0.9 0.3 0.1]); hold on;
else
    plot([0,maxlab],[p(2),(p(2) + maxlab*p(1))],'Color',[0.9 0.3 0.1]); hold on;
end
xlim([0 max(cat(1,bigsst,datsst))]);
ylim([0 max(cat(1,bigsst,datsst))]);

subplot(1,3,3); hold on; box on;
title('Vip','FontSize',20)
plot([1 max(cat(1,bigvip,datvip))],[1 max(cat(1,bigvip,datvip))],'k-'); hold on;
txt = {sprintf('R_C = %.2f',LinRcalc(bigvip,datvip));...
    sprintf('R = %.2f',corr(bigvip,datvip))};
text(2.5*10^4,5.7*10^4,txt,'FontSize',15);
xlabel('Inferred Counts','FontSize',20);
maxlab = max(cat(1,bigvip,datvip));
set(gca,'YTick',[0 0.5*maxlab maxlab]);
set(gca,'YTickLabels',{'0',sprintf('%.1f',0.5*maxlab/10^4),'x10^4'});
set(gca,'XTick',[0 0.5*maxlab maxlab]);
set(gca,'XTickLabels',{'0',sprintf('%.1f',0.5*maxlab/10^4),'x10^4'});
set(gca,'FontSize',18);
p = polyfit(bigvip,datvip,1);
x_maxy = (maxlab - p(2))/p(1);
if x_maxy < maxlab
    plot([0,x_maxy],[p(2),(p(2) + x_maxy*p(1))],'Color',[0.9 0.3 0.1]); hold on;
else
    plot([0,maxlab],[p(2),(p(2) + maxlab*p(1))],'Color',[0.9 0.3 0.1]); hold on;
end
xlim([0 max(cat(1,bigvip,datvip))]);
ylim([0 max(cat(1,bigvip,datvip))]);

if savenclose
    print('Figure_4b_interneuron_scatter','-dtiffn');
    close
end
end