function S_Figure_crossvalcomp(savenclose, directory)

if nargin < 2
    directory = [cd filesep 'MatFiles']; 
    if nargin < 1
        savenclose = 0;
    end
end

LinRcalc = @(x,y) 2*corr(x,y)*std(x)*std(y)/(std(x)^2 + std(y)^2 + (mean(x) - mean(y))^2);
loadstruct = load([directory filesep 'tasic_MRx3_l250_nocv.mat'],'fitstruct');
ngs_nocv = zeros(1,length(loadstruct.fitstruct)); negloglike_nocv = ngs_nocv;
for i = 1:length(ngs_nocv)
    ngs_nocv(i) = loadstruct.fitstruct(i).nGen;
    negloglike_nocv(i) = loadstruct.fitstruct(i).negloglikelihood;
end

loadstruct = load([directory filesep 'tasic_MRx3_l250.mat'],'fitstruct');
ngs_cv = zeros(1,length(loadstruct.fitstruct)); negloglike_cv = ngs_cv;
for i = 1:length(ngs_cv)
    ngs_cv(i) = loadstruct.fitstruct(i).nGen;
    negloglike_cv(i) = loadstruct.fitstruct(i).negloglikelihood;
end

common_ngs = intersect(ngs_cv,ngs_nocv);
negloglike_cv = negloglike_cv(ismember(ngs_cv,common_ngs));
negloglike_nocv = negloglike_nocv(ismember(ngs_nocv,common_ngs));
plotmin = min([negloglike_cv, negloglike_nocv]);
plotmax = max([negloglike_cv, negloglike_nocv]);

figure; hold on;
scatter(negloglike_cv,negloglike_nocv,'bo','filled'); hold on;
plot([plotmin,plotmax],[plotmin,plotmax],'k--','LineWidth',2);
ylim([1.01*plotmin, 0.99*plotmax]); ytickformat('%.2f');
xlim([1.01*plotmin, 0.99*plotmax]); xtickformat('%.2f');
set(gca,'FontSize',18);
set(gca,'XTick',[plotmin (plotmin + plotmax)/2, plotmax]);
set(gca,'YTick',[plotmin (plotmin + plotmax)/2, plotmax]);
legend({'-log{\it p}_{GMM}({\it C}_{red}|{\it n}_G,\lambda)' ,sprintf('R_c = %.2f',LinRcalc(negloglike_cv.',negloglike_nocv.'))},'Location','northwest');
ylabel('No CV','FontSize',20);
xlabel('5-Fold CV','FontSize',20);
title('GMM Fitting With and Without CV','FontSize',20);

if savenclose
    print('sfig_crossvalcomp','-dtiffn');
    close;
end

end