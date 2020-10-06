function Figure_3ab_methodcomp(savenclose, directory)

if nargin < 2
    directory = [cd filesep 'MatFiles']; 
    if nargin < 1
        savenclose = 0;
    end
end
lambda = 250;
methodnames = {'colAMD', 'DBSCAN', 'Entropy', sprintf('MRx3_l%d_coarsesigmas',lambda)};
ngcell = cell(1,length(methodnames)); negloglikescell = ngcell; 
sumfitcell = ngcell; minngcell = ngcell;
for j = 1:length(methodnames)
    loadstruct = load([directory filesep 'tasic_' methodnames{j} '.mat'],'fitstruct','sumfits','minngs');
    ng_list = zeros(1,length(loadstruct.fitstruct));
    negloglike = ng_list;
    for i = 1:length(loadstruct.fitstruct)
        ng_list(i) = loadstruct.fitstruct(i).nGen;
        negloglike(i) = loadstruct.fitstruct(i).negloglikelihood;
    end
    ngcell{j} = ng_list;
    negloglikescell{j} = negloglike;
    sumfitcell{j} = loadstruct.sumfits;
    minngcell{j} = loadstruct.minngs;
end

minnegloglikes = zeros(1,length(negloglikescell));
maxnegloglikes = minnegloglikes;
maxsf = minnegloglikes; minsf = maxsf;
for i = 1:length(negloglikescell)
    minnegloglikes(i) = min(negloglikescell{i});
    maxnegloglikes(i) = max(negloglikescell{i});
    maxsf(i) = max(sumfitcell{i});
    minsf(i) = min(sumfitcell{i});
end

f2 = figure; hold on;
set(f2,'Position',[0 0 500 1000]); hold on;
figure(f2);
subplot(2,1,1)
for i = 1:length(methodnames)
    plot(ngcell{i},negloglikescell{i},'LineWidth',2.5); hold on;
end
legend({'colAMD', 'DBSCAN', 'Entropy', sprintf('MRx3, lambda = %d',lambda)});
ylim([1.1*min(minnegloglikes), 0.9*max(maxnegloglikes)]); ytickformat('%.2f');
xlim([0 3855]);
set(gca,'FontSize',18);
set(gca,'XTick',[0 1750 3500]);
set(gca,'YTick',[min(minnegloglikes), (min(minnegloglikes)+max(maxnegloglikes))/2,...
    max(maxnegloglikes)]);
ylabel('-log{\it p}_{GMM}({\it C}_{red}|{\it n}_G,\lambda)','FontSize',20);
xlabel('{\it n}_G','FontSize',20);
title('Method Comparison, Classification','FontSize',20);
hold off;
subplot(2,1,2);
for i = 1:length(methodnames)
    plot(minngcell{i},sumfitcell{i},'LineWidth',2.5); hold on;
end
ylim([0.8*min(minsf), 1.2*max(maxsf)]); ytickformat('%.2f');
xlim([400 850]);
legend({'colAMD', 'DBSCAN', 'Entropy', sprintf('MRx3, lambda = %d',lambda)});
set(gca,'FontSize',18);
% set(gca,'XTick',[3]);
set(gca,'YTick',[min(minsf), (min(minsf)+max(maxsf))/2,  max(maxsf)]);
ylabel('\Sigma_{fit}','FontSize',20);
xlabel('{\it n}_G^*','FontSize',20);
% title({'Method','Comparison'},'FontSize',25);
title('Method Comparison, \Sigma_{fit}','FontSize',20);

if savenclose
    print('s_figure_methodcomp','-dtiffn')
    close
end

end