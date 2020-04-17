function BlandAltman_figplot_micro(outstruct,idx)
% Designed to handle an outstruct of dimension length(ngenelist) that has a
% single field "bvals." idx is the index of interest that is within the
% range of ngenelist
ngenes = outstruct(idx).nGen;
neoinds = 57:94; 
foreinds = [1:11,23:56,141:156,170:212];
wbinds = [12:22 95:140 157:169];
indcell = {wbinds,foreinds,neoinds};
indcolor = {'b','m','g'};
% indmark = {'x','+','d','*','o','s'};
indmark = {'*','o','d','s','^','p'};
% reslabels = {'Hind & Midbrain','Other Forebrain','Neocortex'};
celllabs = {'Pv','Sst','Vip','Total','Micro.','Neuron'};
LinRcalc = @(x,y) 2*corr(x,y)*std(x)*std(y)/(std(x)^2 + std(y)^2 + (mean(x) - mean(y))^2);
bigmicro = 0; datmicro = 0;

f1 = figure('Units','inch','Position',[0 0 16 4]); hold on;

for i = 1:length(indcell)
    load('keller_micro_totals_reorder.mat');
    load('keller_micro_listB.mat');
    
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
    
    cellsim = infer_micro.';
    celldat = hh_micro.';
    
    bigmicro = [bigmicro; infer_micro.'];
    datmicro = [datmicro; hh_micro.'];
    
    set(gca,'FontSize',18);
    scatter(cellsim_cell{5},celldat_cell{5},60,[indcolor{i} indmark{5}],'filled'); hold on;
    ylabel([celllabs{5} ' Counts'],'FontSize',22)
    xlabel('Inferred Counts','FontSize',22)
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

% subplot(15,5,[36 37 41 42]); hold on; box on;
% title([sprintf('R_C = %.2f',LinRcalc(bigpv,datpv)) ', ' sprintf('R = %.2f',corr(bigpv,datpv))],'FontSize',12)
% plot([1 max([bigpv;datpv])],[1 max(cat(1,bigpv,datpv))],'k-'); hold on;
% % txt = {sprintf('R_C = %.2f',LinRcalc(bigpv,datpv));...
% %     sprintf('R = %.2f',corr(bigpv,datpv))};
% % if strcmp(method,'mRMR')
% %     text(0.1*10^4,5.5*10^4,txt,'FontSize',8);
% % elseif strcmp(method,'SoRDES')
% %     text(0.1*10^4,7.5*10^4,txt,'FontSize',8);
% % end
% maxlab = max(cat(1,bigpv,datpv));
% set(gca,'YTick',[0 0.5*maxlab maxlab]);
% set(gca,'YTickLabels',{'0',sprintf('%.1f',0.5*maxlab/10^5),'x10^5'});
% set(gca,'XTick',[0 0.5*maxlab maxlab]);
% set(gca,'XTickLabels',{'0',sprintf('%.1f',0.5*maxlab/10^5),'x10^5'});
% p = polyfit(bigpv,datpv,1);
% % y = polyval(p,bigpv);
% x_maxy = (maxlab - p(2))/p(1);
% if x_maxy < maxlab
%     plot([0,x_maxy],[p(2),(p(2) + x_maxy*p(1))],'Color',[0.9 0.3 0.1]); hold on;
% else
%     plot([0,maxlab],[p(2),(p(2) + maxlab*p(1))],'Color',[0.9 0.3 0.1]); hold on;
% end
% % plot(bigpv,y,'Color',[0.9 0.3 0.1]); hold on;
% xlim([0 max(cat(1,bigpv,datpv))]);
% ylim([0 max(cat(1,bigpv,datpv))]);
% 
% subplot(15,5,[39 40 44 45]); hold on; hold on; box on;
% title([sprintf('R_C = %.2f',LinRcalc(bigsst,datsst)) ', ' sprintf('R = %.2f',corr(bigsst,datsst))],'FontSize',12)
% plot([1 max(cat(1,bigsst,datsst))],[1 max(cat(1,bigsst,datsst))],'k-'); hold on;
% % txt = {sprintf('R_C = %.2f',LinRcalc(bigsst,datsst));...
% %     sprintf('R = %.2f',corr(bigsst,datsst))};
% % text(3.5*10^4,1*10^4,txt,'FontSize',8);
% maxlab = max(cat(1,bigsst,datsst));
% set(gca,'YTick',[0 0.5*maxlab maxlab]);
% set(gca,'YTickLabels',{'0',sprintf('%.1f',0.5*maxlab/10^4),'x10^4'});
% set(gca,'XTick',[0 0.5*maxlab maxlab]);
% set(gca,'XTickLabels',{'0',sprintf('%.1f',0.5*maxlab/10^4),'x10^4'});
% p = polyfit(bigsst,datsst,1);
% x_maxy = (maxlab - p(2))/p(1);
% if x_maxy < maxlab
%     plot([0,x_maxy],[p(2),(p(2) + x_maxy*p(1))],'Color',[0.9 0.3 0.1]); hold on;
% else
%     plot([0,maxlab],[p(2),(p(2) + maxlab*p(1))],'Color',[0.9 0.3 0.1]); hold on;
% end
% % y = polyval(p,bigsst);
% % plot(bigsst,y,'Color',[0.9 0.3 0.1]); hold on;
% xlim([0 max(cat(1,bigsst,datsst))]);
% ylim([0 max(cat(1,bigsst,datsst))]);
% 
% subplot(15,5,[51 52 56 57]); hold on; box on;
% title([sprintf('R_C = %.2f',LinRcalc(bigvip,datvip)) ', ' sprintf('R = %.2f',corr(bigvip,datvip))],'FontSize',12)
% plot([1 max(cat(1,bigvip,datvip))],[1 max(cat(1,bigvip,datvip))],'k-'); hold on;
% % txt = {sprintf('R_C = %.2f',LinRcalc(bigvip,datvip));...
% %     sprintf('R = %.2f',corr(bigvip,datvip))};
% % if strcmp(method,'mRMR')
% %     text(0.1*10^4,3.5*10^4,txt,'FontSize',8);
% % elseif strcmp(method,'SoRDES')
% %     text(0.1*10^4,2.25*10^4,txt,'FontSize',8);
% % end
% maxlab = max(cat(1,bigvip,datvip));
% set(gca,'YTick',[0 0.5*maxlab maxlab]);
% set(gca,'YTickLabels',{'0',sprintf('%.1f',0.5*maxlab/10^4),'x10^4'});
% set(gca,'XTick',[0 0.5*maxlab maxlab]);
% set(gca,'XTickLabels',{'0',sprintf('%.1f',0.5*maxlab/10^4),'x10^4'});
% p = polyfit(bigvip,datvip,1);
% x_maxy = (maxlab - p(2))/p(1);
% if x_maxy < maxlab
%     plot([0,x_maxy],[p(2),(p(2) + x_maxy*p(1))],'Color',[0.9 0.3 0.1]); hold on;
% else
%     plot([0,maxlab],[p(2),(p(2) + maxlab*p(1))],'Color',[0.9 0.3 0.1]); hold on;
% end
% % y = polyval(p,bigvip);
% % plot(bigvip,y,'Color',[0.9 0.3 0.1]); hold on;
% xlim([0 max(cat(1,bigvip,datvip))]);
% ylim([0 max(cat(1,bigvip,datvip))]);
% 
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

% subplot(15,5,[66 67 71 72]); 
hold on; box on;
title('Microglia','FontSize',22)
plot([1 max(cat(1,bigmicro,datmicro))],[1 max(cat(1,bigmicro,datmicro))],'k-'); hold on;
txt = {sprintf('R_C = %.2f',LinRcalc(bigmicro,datmicro));...
    sprintf('R = %.2f',corr(bigmicro,datmicro))};
text(2.2*10^5,0.5*10^5,txt,'FontSize',18);
% txt = {sprintf('R_C = %.2f',LinRcalc(bigmicro,datmicro));...
%     sprintf('R = %.2f',corr(bigmicro,datmicro))};
% if strcmp(method,'mRMR')
%     text(0.1*10^5,4.75*10^5,txt,'FontSize',8);
% elseif strcmp(method,'SoRDES')
%     text(0.1*10^5,2.5*10^5,txt,'FontSize',8);
% end
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
% y = polyval(p,bigmicro);
% plot(bigmicro,y,'Color',[0.9 0.3 0.1]); hold on;
xlim([0 max(cat(1,bigmicro,datmicro))]);
ylim([0 max(cat(1,bigmicro,datmicro))]);

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