function S_Figure_MRx3_acrossnG(outstruct,nG_opt,savenclose)

if nargin < 3
    savenclose = 0;
    if nargin < 2
        nG_opt = 529;
    end
end

for i = 1:length(outstruct)
    ngs(i) = outstruct(i).nGen;
    sumfits(i) = outstruct(i).sumfit;
    pvrc(i) = outstruct(i).LinR.pv(end);
    viprc(i) = outstruct(i).LinR.vip(end);
    sstrc(i) = outstruct(i).LinR.sst(end);
    microrc(i) = outstruct(i).LinR.micro(end);
    taus(i) = outstruct(i).tau;
end

cmap = lines(6);
f2 = figure;
set(f2,'Position',[0 0 600 550]);
sgtitle({'\bf MRx3, \lambda = 250, Tasic et al. Comparisons'},'FontSize', 25);
subplot(3,3,[1,2,4,5])
hold on;
plot(ngs,sumfits,'Color',cmap(1,:),'LineWidth',3);
plot([nG_opt nG_opt], [0.9*min(sumfits) 1.1*max(sumfits)],'k--','LineWidth',2);
xlim([0 3855]); xticks([0 1750 3500]);
ylim([0.9*min(sumfits) 1.1*max(sumfits)]); yticks([min(sumfits) mean([max(sumfits),min(sumfits)]) max(sumfits)]);
ytickformat('%.2f'); set(gca,'FontSize',14);
xlabel('{\it n}_G','FontSize',18); ylabel('\Sigma_{fit}','FontSize',18);
hold off;

subplot(3,3,3)
hold on;
plot(ngs,pvrc,'Color',cmap(2,:),'LineWidth',2.5);
plot([nG_opt nG_opt], [0.9*min(pvrc) 1.1*max(pvrc)],'k--','LineWidth',1.5);
xlim([0 3855]); xticks([0 1750 3500]);
ylim([0.9*min(pvrc) 1.1*max(pvrc)]); yticks([min(pvrc) mean([max(pvrc),min(pvrc)]) max(pvrc)]);
ytickformat('%.2f'); set(gca,'FontSize',14);
xlabel('{\it n}_G','FontSize',16); ylabel('{\it Pvalb}+ R_c','FontSize',16);
hold off;

subplot(3,3,6)
hold on;
plot(ngs,sstrc,'Color',cmap(3,:),'LineWidth',2.5);
plot([nG_opt nG_opt], [0.9*min(sstrc) 1.1*max(sstrc)],'k--','LineWidth',1.5);
xlim([0 3855]); xticks([0 1750 3500]);
ylim([0.9*min(sstrc) 1.1*max(sstrc)]); yticks([min(sstrc) mean([max(sstrc),min(sstrc)]) max(sstrc)]);
ytickformat('%.2f'); set(gca,'FontSize',14);
xlabel('{\it n}_G','FontSize',16); ylabel('{\it Sst}+ R_c','FontSize',16);
hold off;

subplot(3,3,7)
hold on;
plot(ngs,taus,'Color',cmap(4,:),'LineWidth',2.5);
plot([nG_opt nG_opt], [0.9*min(taus) 1.1*max(taus)],'k--','LineWidth',1.5);
xlim([0 3855]); xticks([0 1750 3500]);
ylim([0.9*min(taus) 1.1*max(taus)]); yticks([min(taus) mean([max(taus),min(taus)]) max(taus)]);
ytickformat('%.2f'); set(gca,'FontSize',14);
xlabel('{\it n}_G','FontSize',16); ylabel('\tau_{adj}','FontSize',16);
hold off;

subplot(3,3,8)
hold on;
plot(ngs,microrc,'Color',cmap(5,:),'LineWidth',2.5);
plot([nG_opt nG_opt], [0.9*min(microrc) 1.1*max(microrc)],'k--','LineWidth',1.5);
xlim([0 3855]); xticks([0 1750 3500]);
ylim([0.9*min(microrc) 1.1*max(microrc)]); yticks([min(microrc) mean([max(microrc),min(microrc)]) max(microrc)]);
ytickformat('%.2f'); set(gca,'FontSize',14);
xlabel('{\it n}_G','FontSize',16); ylabel('Micro R_c','FontSize',16);
hold off;

subplot(3,3,9)
hold on;
plot(ngs,viprc,'Color',cmap(6,:),'LineWidth',2.5);
plot([nG_opt nG_opt], [0.9*min(viprc) 1.1*max(viprc)],'k--','LineWidth',1.5);
xlim([0 3855]); xticks([0 1750 3500]);
ylim([0.9*min(viprc) 1.1*max(viprc)]); yticks([min(viprc) mean([max(viprc),min(viprc)]) max(viprc)]);
ytickformat('%.2f'); set(gca,'FontSize',14);
xlabel('{\it n}_G','FontSize',16); ylabel('{\it Vip}+ R_c','FontSize',16);
hold off;

if savenclose
    print('mrx3acrossng','-dtiffn');
    close;
end
end