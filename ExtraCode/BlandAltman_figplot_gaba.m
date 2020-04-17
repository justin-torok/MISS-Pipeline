function [allcellsmat,bigsim,bigdat] = Figure_3b_interneuron(outstruct,idx)
% Designed to handle an outstruct of dimension length(ngenelist) that has a
% single field "bvals." idx is the index of interest that is within the
% range of ngenelist
ngenes = outstruct(idx).nGen;
sampleinds = [70,92,186]; % Secondary motor area, primary visual area, dorsal lateral geniculate 
neoinds = [57:69, 71:91, 93, 94]; 
foreinds = [1:11,23:56,141:156,170:185, 187:212];
wbinds = [12:22 95:140 157:169];
indcell = {wbinds,foreinds,neoinds,sampleinds};
indcolor = {'b','m','g','r'};
indsize =  [30,30,30,100];
% indmark = {'x','+','d','*','o','s'};
indmark = {'s','o','d','*','^','p'};
reslabels = {'Hind & Midbrain','Other Forebrain','Neocortex','Sampled Regions'};
celllabs = {'Pv','Sst','Vip','Total','Micro.','Neuron'};

LinRcalc = @(x,y) 2*corr(x,y)*std(x)*std(y)/(std(x)^2 + std(y)^2 + (mean(x) - mean(y))^2);

% BA_diffs_kim = cell(3,length(indcell));
% BA_avgs_kim = BA_diffs_kim;
% BA_diffs_murakami = cell(length(indcell));
% BA_avgs_murakami = BA_diffs_murakami;

bigsim = 0; bigdat = 0;
bigpv = 0; datpv = 0;
bigsst = 0; datsst = 0;
bigvip = 0; datvip = 0;
bigall = 0; datall = 0;
bigmicro = 0; datmicro = 0;
bigneuron = 0; datneuron = 0;
f1 = figure('Units','inch','Position',[0 0 16 4]);
% set(f1,'Position',[0 0 400 1000]); hold on; %this to set the size
s1 = subplot(3,1,1,'Parent',f1); hold on; %this to make it a subplot at a certain size
for i = 1:length(indcell)
    load('kim_totals_reorder_m.mat');
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

    load('murakami_totals_reorder.mat');
    murakami_totals_reorder = murakami_totals_reorder(testinds);
    nonzeroinds = find(murakami_totals_reorder);
    murakami_tot = murakami_totals_reorder(nonzeroinds);

    sumcts_wb = outstruct(idx).Bsums(1:213,:) + outstruct(idx).Bsums(214:end,:);
    sumcts_wb(12,:) = [];
    sumcts_wb_mod = sumcts_wb(testinds,:);
    sumcts_pv_wb = sumcts_wb_mod(nonnaninds(:,1),3);
    sumcts_sst_wb = sumcts_wb_mod(nonnaninds(:,2),6);
    sumcts_vip_wb = sumcts_wb_mod(nonnaninds(:,3),7);
    sumcts_sum = sum(sumcts_wb_mod(nonzeroinds,:),2);
    
    load('herculano_houzel_micro_totals_reorder.mat');
    load('herculano_houzel_micro_listB.mat');
    load('herculano_houzel_neuron_totals_reorder.mat');
    load('herculano_houzel_neuron_listB.mat');
    
    sumcts_micro_wb = sumcts_wb(:,23);
    sumcts_neuron_wb = sum(sumcts_wb(:,1:21),2);
    
    hh_micro = [];
    infer_micro = [];
    for j = 1:size(herculano_houzel_micro_listB)
        hh = herculano_houzel_micro_listB(j,:);
        hhinds = hh(hh > 0);
        if all(ismember(hhinds,testinds))
            hh_micro = [hh_micro, sum(hh_micro_totals_reorder(hhinds))];
            infer_micro = [infer_micro, sum(sumcts_micro_wb(hhinds))];
        end
    end
    
    hh_neuron = [];
    infer_neuron = [];
    for j = 1:size(herculano_houzel_neuron_listB)
        hh = herculano_houzel_neuron_listB(j,:);
        hhinds = hh(hh > 0);
        if all(ismember(hhinds,testinds))
            hh_neuron = [hh_neuron, sum(hh_neuron_totals_reorder(hhinds))];
            infer_neuron = [infer_neuron, sum(sumcts_neuron_wb(hhinds))];
        end
    end
    
    if isempty(infer_neuron)
        infer_neuron = 0;
    end
    if isempty(hh_neuron)
        hh_neuron = 0;
    end
    
    cellsim = [sumcts_pv_wb; sumcts_sst_wb; sumcts_vip_wb;...
        sumcts_sum; infer_micro.'; infer_neuron.'];
    
    celldat = [kim_pv_wb; kim_sst_wb; kim_vip_wb;...
        murakami_tot; hh_micro.'; hh_neuron.'];
    
    cellsim_cell = {sumcts_pv_wb, sumcts_sst_wb, sumcts_vip_wb,...
                    sumcts_sum, infer_micro.', infer_neuron.'};
                
    celldat_cell = {kim_pv_wb, kim_sst_wb, kim_vip_wb,...
                    murakami_tot, hh_micro.', hh_neuron.'};
                
    allcellsmat(i).sim = cellsim;
    allcellsmat(i).dat = celldat;
    
    bigdat = [bigdat; celldat];
    bigsim = [bigsim; cellsim];
    bigpv = [bigpv; sumcts_pv_wb];
    datpv = [datpv; kim_pv_wb];
    bigsst = [bigsst; sumcts_sst_wb];
    datsst = [datsst; kim_sst_wb];
    bigvip = [bigvip; sumcts_vip_wb];
    datvip = [datvip; kim_vip_wb];
    bigall = [bigall; sumcts_sum];
    datall = [datall; murakami_tot];
    bigmicro = [bigmicro; infer_micro.'];
    datmicro = [datmicro; hh_micro.'];
    bigneuron = [bigneuron; infer_neuron.'];
    datneuron = [datneuron; hh_neuron.'];
    
    allcellsmat(i).LinR = LinRcalc(cellsim,celldat);
    allcellsmat(i).RVal = corr(cellsim,celldat);
    
%     for j = 1:length(indmark)
%         if ismember(j,[1])
%             subplot(15,5,1:25); hold on;
%             set(gca,'FontSize',14);
%             scatter(cellsim_cell{j},celldat_cell{j},50,[indcolor{i} indmark{j}]);
%         else
%             subplot(15,5,1:25); hold on;
%             set(gca,'FontSize',14);
%             scatter(cellsim_cell{j},celldat_cell{j},50,[indcolor{i} indmark{j}],'filled');
%         end
%     end
    
    subplot(1,3,1); hold on;
    set(gca,'FontSize',12);
    scatter(cellsim_cell{1},celldat_cell{1},indsize(i),[indcolor{i} indmark{1}],'filled'); hold on;
%     ylabel([celllabs{1} ' Counts'],'FontSize',14)
    
    subplot(1,3,2); hold on;
    set(gca,'FontSize',12);
    scatter(cellsim_cell{2},celldat_cell{2},indsize(i),[indcolor{i} indmark{2}],'filled'); hold on;
%     ylabel([celllabs{2} ' Counts'],'FontSize',14)
    
    subplot(1,3,3); hold on;
    set(gca,'FontSize',12);
    scatter(cellsim_cell{3},celldat_cell{3},indsize(i),[indcolor{i} indmark{3}],'filled'); hold on;
%     ylabel([celllabs{3} ' Counts'],'FontSize',14)
    
%     subplot(15,5,[54 55 59 60]); hold on;
%     set(gca,'FontSize',12);
%     scatter(cellsim_cell{4},celldat_cell{4},30,[indcolor{i} indmark{4}],'filled'); hold on;
%     ylabel([celllabs{4} ' Counts'],'FontSize',14) 
%     
%     subplot(15,5,[66 67 71 72]); hold on;
%     set(gca,'FontSize',12);
%     scatter(cellsim_cell{5},celldat_cell{5},30,[indcolor{i} indmark{5}],'filled'); hold on;
%     ylabel([celllabs{5} ' Counts'],'FontSize',14)
%     xlabel('Inferred Counts','FontSize',14)
%     
%     subplot(15,5,[69 70 74 75]); hold on;
%     set(gca,'FontSize',12);
%     scatter(cellsim_cell{6},celldat_cell{6},30,[indcolor{i} indmark{6}],'filled'); hold on;
%     ylabel([celllabs{6} ' Counts'],'FontSize',14)
%     xlabel('Inferred Counts','FontSize',14)

end

bigpv(1) = []; bigsst(1) = []; bigvip(1) = []; bigall(1) = [];
bigmicro(1) = []; bigneuron(1) = [];
datpv(1) = []; datsst(1) = []; datvip(1) = []; datall(1) = [];
datmicro(1) = []; datneuron(1) = [];

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
% y = polyval(p,bigpv);
x_maxy = (maxlab - p(2))/p(1);
if x_maxy < maxlab
    plot([0,x_maxy],[p(2),(p(2) + x_maxy*p(1))],'Color',[0.9 0.3 0.1]); hold on;
else
    plot([0,maxlab],[p(2),(p(2) + maxlab*p(1))],'Color',[0.9 0.3 0.1]); hold on;
end
% plot(bigpv,y,'Color',[0.9 0.3 0.1]); hold on;
xlim([0 max(cat(1,bigpv,datpv))]);
ylim([0 max(cat(1,bigpv,datpv))]);

subplot(1,3,2); hold on; hold on; box on;
title('Sst','FontSize',20)
plot([1 max(cat(1,bigsst,datsst))],[1 max(cat(1,bigsst,datsst))],'k-'); hold on;
txt = {sprintf('R_C = %.2f',LinRcalc(bigsst,datsst));...
    sprintf('R = %.2f',corr(bigsst,datsst))};
text(3.75*10^4,8.5*10^4,txt,'FontSize',15);
% ylabel('Cell Counts','FontSize',20);
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
% y = polyval(p,bigsst);
% plot(bigsst,y,'Color',[0.9 0.3 0.1]); hold on;
xlim([0 max(cat(1,bigsst,datsst))]);
ylim([0 max(cat(1,bigsst,datsst))]);

subplot(1,3,3); hold on; box on;
title('Vip','FontSize',20)
plot([1 max(cat(1,bigvip,datvip))],[1 max(cat(1,bigvip,datvip))],'k-'); hold on;
txt = {sprintf('R_C = %.2f',LinRcalc(bigvip,datvip));...
    sprintf('R = %.2f',corr(bigvip,datvip))};
text(2.5*10^4,5.7*10^4,txt,'FontSize',15);
% ylabel('Cell Counts','FontSize',20)
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
% y = polyval(p,bigvip);
% plot(bigvip,y,'Color',[0.9 0.3 0.1]); hold on;
xlim([0 max(cat(1,bigvip,datvip))]);
ylim([0 max(cat(1,bigvip,datvip))]);

% subplot(15,5,[54 55 59 60]); hold on; box on;
% title([sprintf('R_C = %.2f',LinRcalc(bigall,datall)) ', ' sprintf('R = %.2f',corr(bigall,datall))],'FontSize',12)
% plot([1 max(cat(1,bigall,datall))],[1 max(cat(1,bigall,datall))],'k-'); hold on;
% % txt = {sprintf('R_C = %.2f',LinRcalc(bigall,datall));...
% %     sprintf('R = %.2f',corr(bigall,datall))};
% % text(10*10^5,3*10^5,txt,'FontSize',8);
% maxlab = max(cat(1,bigall,datall));
% set(gca,'YTick',[0 0.5*maxlab maxlab]);
% set(gca,'YTickLabels',{'0',sprintf('%.1f',0.5*maxlab/10^6),'x10^6'});
% set(gca,'XTick',[0 0.5*maxlab maxlab]);
% set(gca,'XTickLabels',{'0',sprintf('%.1f',0.5*maxlab/10^6),'x10^6'});
% p = polyfit(bigall,datall,1);
% x_maxy = (maxlab - p(2))/p(1);
% if x_maxy < maxlab
%     plot([0,x_maxy],[p(2),(p(2) + x_maxy*p(1))],'Color',[0.9 0.3 0.1]); hold on;
% else
%     plot([0,maxlab],[p(2),(p(2) + maxlab*p(1))],'Color',[0.9 0.3 0.1]); hold on;
% end
% % y = polyval(p,bigall);
% % plot(bigall,y,'Color',[0.9 0.3 0.1]); hold on;
% xlim([0 max(cat(1,bigall,datall))]);
% ylim([0 max(cat(1,bigall,datall))]);
% 
% subplot(15,5,[66 67 71 72]); hold on; box on;
% title([sprintf('R_C = %.2f',LinRcalc(bigmicro,datmicro)) ', ' sprintf('R = %.2f',corr(bigmicro,datmicro))],'FontSize',12)
% plot([1 max(cat(1,bigmicro,datmicro))],[1 max(cat(1,bigmicro,datmicro))],'k-'); hold on;
% % txt = {sprintf('R_C = %.2f',LinRcalc(bigmicro,datmicro));...
% %     sprintf('R = %.2f',corr(bigmicro,datmicro))};
% % if strcmp(method,'mRMR')
% %     text(0.1*10^5,4.75*10^5,txt,'FontSize',8);
% % elseif strcmp(method,'SoRDES')
% %     text(0.1*10^5,2.5*10^5,txt,'FontSize',8);
% % end
% maxlab = max(cat(1,bigmicro,datmicro));
% set(gca,'YTick',[0 0.5*maxlab maxlab]);
% set(gca,'YTickLabels',{'0',sprintf('%.1f',0.5*maxlab/10^5),'x10^5'});
% set(gca,'XTick',[0 0.5*maxlab maxlab]);
% set(gca,'XTickLabels',{'0',sprintf('%.1f',0.5*maxlab/10^5),'x10^5'});
% p = polyfit(bigmicro,datmicro,1);
% x_maxy = (maxlab - p(2))/p(1);
% if x_maxy < maxlab
%     plot([0,x_maxy],[p(2),(p(2) + x_maxy*p(1))],'Color',[0.9 0.3 0.1]); hold on;
% else
%     plot([0,maxlab],[p(2),(p(2) + maxlab*p(1))],'Color',[0.9 0.3 0.1]); hold on;
% end
% % y = polyval(p,bigmicro);
% % plot(bigmicro,y,'Color',[0.9 0.3 0.1]); hold on;
% xlim([0 max(cat(1,bigmicro,datmicro))]);
% ylim([0 max(cat(1,bigmicro,datmicro))]);
% 
% subplot(15,5,[69 70 74 75]); hold on; box on;
% title([sprintf('R_C = %.2f',LinRcalc(bigneuron,datneuron)) ', ' sprintf('R = %.2f',corr(bigneuron,datneuron))],'FontSize',12)
% plot([1 max(cat(1,bigneuron,datneuron))],[1 max(cat(1,bigneuron,datneuron))],'k-'); hold on;
% % txt = {sprintf('R_C = %.2f',LinRcalc(bigneuron,datneuron));...
% %     sprintf('R = %.2f',corr(bigneuron,datneuron))};
% % if strcmp(method,'mRMR')
% %     text(0.25*10^5,15.5*10^5,txt,'FontSize',8);
% % elseif strcmp(method,'SoRDES')
% %     text(0.25*10^5,15*10^5,txt,'FontSize',8);
% % end
% maxlab = max(cat(1,bigneuron,datneuron));
% set(gca,'YTick',[0 0.5*maxlab maxlab]);
% set(gca,'YTickLabels',{'0',sprintf('%.1f',0.5*maxlab/10^5),'x10^5'});
% set(gca,'XTick',[0 0.5*maxlab maxlab]);
% set(gca,'XTickLabels',{'0',sprintf('%.1f',0.5*maxlab/10^5),'x10^5'});
% p = polyfit(bigneuron,datneuron,1);
% x_maxy = (maxlab - p(2))/p(1);
% if x_maxy < maxlab
%     plot([0,x_maxy],[p(2),(p(2) + x_maxy*p(1))],'Color',[0.9 0.3 0.1]); hold on;
% else
%     plot([0,maxlab],[p(2),(p(2) + maxlab*p(1))],'Color',[0.9 0.3 0.1]); hold on;
% end
% % y = polyval(p,bigneuron);
% % plot(bigneuron,y,'Color',[0.9 0.3 0.1]); hold on;
% xlim([0 max(cat(1,bigneuron,datneuron))]);
% ylim([0 max(cat(1,bigneuron,datneuron))]);
% 
% bigdat(1) = []; bigsim(1) = [];
% subplot(15,5,1:25); hold on; box on;
% plot([1 max(cat(1,bigsim,bigdat))],[1 max(cat(1,bigsim,bigdat))],'k-'); hold on;
% ylim([1 max(bigdat)]); xlim([1 max(bigdat)]);
% yticks([10^0 10^2 10^4 10^6]);
% xticks([10^0 10^2 10^4 10^6]);
% set(gca,'XScale','log');
% set(gca,'Yscale','log');
% set(gca,'FontSize',14);
% ylabel('Literature Counts','FontSize',16);
% xlabel('Inferred Counts','FontSize',16);
% title([method ', Literature vs. Inferred'],'FontSize',20);
% txt = {['{\it n}_G ' sprintf('= %d',ngenes)];...
%     sprintf('R_C = %.2f',LinRcalc(bigdat,bigsim));...
%     sprintf('R = %.2f',corr(bigdat,bigsim))};
% text(0.125*10^1,4.85*10^5,txt,'FontSize',10);

end