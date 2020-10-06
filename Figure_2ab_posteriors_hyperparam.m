function Figure_2ab_posteriors_hyperparam(savenclose, directory)

if nargin < 2
    directory = [cd filesep 'MatFiles']; 
    if nargin < 1
        savenclose = 0;
    end
end

lambda = 250; % hard-coded
sigmas = [3400,4400,6000,8000]; % hard-coded

method = sprintf('MRx3_l%d',lambda);
loadstruct = load([directory filesep 'tasic_' method '.mat'],'fitstruct');
ng_list = zeros(1,length(loadstruct.fitstruct));
negloglike = ng_list;
for i = 1:length(loadstruct.fitstruct)
    ng_list(i) = loadstruct.fitstruct(i).nGen;
    negloglike(i) = loadstruct.fitstruct(i).negloglikelihood;
end
neglogpostsmat = zeros(length(sigmas),length(ng_list));
for k = 1:length(sigmas)
    neglogprior = (ng_list.^2)/(2*sigmas(k)^2);
    neglogpostsmat(k,:) = negloglike + neglogprior;
end
minneglogposts = zeros(1,length(sigmas));
maxneglogposts = minneglogposts; bestngs = minneglogposts;
for k = 1:length(sigmas)
    [minneglogposts(k),minind] = min(neglogpostsmat(k,:));
    maxneglogposts(k) = max(neglogpostsmat(k,:));
    bestngs(k) = ng_list(minind);
end
cmap = hsv(length(sigmas));
leg = cell(1,length(sigmas));

% For Figure 2c
lambdas = 0:50:350;
sigmas_ = 3000:25:6000;
sumfitmat = zeros(length(lambdas),length(sigmas_));
for j = 1:length(lambdas)
    load([directory filesep sprintf('tasic_MRx3_l%d.mat',lambdas(j))],'sumfits');
    sumfitmat(j,:) = sumfits;
end

f2 = figure; hold on;
set(f2,'Position',[0 0 500 450]); hold on;
figure(f2);
% subplot(20,1,1:7)
for k = 1:length(sigmas)
    if sigmas(k) == 4400
        lw = 3;
    else
        lw = 1.25;
    end
    plot(ng_list,neglogpostsmat(k,:),'Color', cmap(k,:), 'LineWidth', lw);
end
for k = 1:length(sigmas)
    if sigmas(k) == 4400
        lw = 3;
        ls = '-';
    else
        lw = 1.25;
        ls = '--';
    end
    plot([bestngs(k) bestngs(k)],[1.1*min(minneglogposts), 0.9*max(maxneglogposts)],...
        'Color',cmap(k,:),'LineWidth',lw-0.5,'LineStyle',ls);
    leg{k} = [sprintf('\\sigma = %d, ',sigmas(k)) sprintf('n_G^* = %d',bestngs(k))];
end
ylim([1.1*min(minneglogposts), 0.9*max(maxneglogposts)]); ytickformat('%.2f');
xlim([0 3855]);
legend(leg,'Location','north');
% ylim([0 3.2])
set(gca,'FontSize',20);
set(gca,'XTick',[0 1750 3500]);
set(gca,'YTick',[min(minneglogposts), (min(minneglogposts)+max(maxneglogposts))/2,...
    max(maxneglogposts)]);
ylabel('-log{\it L}({\it n}_G|\lambda,\sigma)','FontSize',25);
xlabel('{\it n}_G','FontSize',25);
% title({'Method','Comparison'},'FontSize',25);
title('Posterior Minimization','FontSize',25);
hold off; 
if savenclose
    print('Figure_2b_posteriormins','-dtiffn')
    close
end

f1 = figure; hold on;
set(f1,'Position',[0 0 500 500]);
% subplot(20,1,11:20)
imagesc(flipud(sumfitmat.')); colormap redblue; colorbar;
set(gca,'FontSize',20,'TickLength',[0 0]);
xlabel('\lambda','FontSize',25); ylabel('\sigma','FontSize',25);
xticks(1:length(lambdas)); xticklabels(num2cell(lambdas)); 
yticklabels(num2cell(fliplr(sigmas_(1:40:121)))); yticks(1:40:121);
ylim([0.5,121.5]); xlim([0.5,8.5]);
xtickangle(45);
title({'Hyperparameter','Optimization'},'FontSize',25);
% title('Hyperparameter Optimization','FontSize',20);

if savenclose
    print('Figure_2c_hyperparam','-dtiffn')
    close
end
end