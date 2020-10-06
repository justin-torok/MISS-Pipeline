function Figure_2c_zeiselposterior(savenclose, directory)

if nargin < 2
    directory = [cd filesep 'MatFiles']; 
    if nargin < 1
        savenclose = 0;
    end
end

lambda = 250; % hard-coded
sigma = 4400;

initlen = 0;
for j = 1:4
    loadstruct = load([directory filesep sprintf('zeisel_gmm_l250_%d.mat',j)],'fitstruct');
    for i = 1:length(loadstruct.fitstruct)
        ng_list(i+initlen) = loadstruct.fitstruct(i).nGen;
        negloglike(i+initlen) = loadstruct.fitstruct(i).negloglikelihood;
    end
    initlen = length(ng_list);
end

neglogprior = (ng_list.^2)/(2*sigma^2);
neglogposts = negloglike + neglogprior;

[minneglogposts,minind] = min(neglogposts);
maxneglogposts = max(neglogposts);
bestngs = ng_list(minind);

leg = {[sprintf('\\sigma = %d, ',sigma) sprintf('n_G^* = %d',bestngs)]};

f2 = figure; hold on;
set(f2,'Position',[0 0 500 450]); hold on;
figure(f2);
lw = 3;
ls = '--';
plot(ng_list,neglogposts,'Color', [1 0 0], 'LineWidth', lw);
plot([bestngs bestngs],[1.1*minneglogposts, 0.9*maxneglogposts],...
    'Color',[0 0 0],'LineWidth',lw-0.5,'LineStyle',ls);
ylim([1.1*minneglogposts, 0.9*maxneglogposts]); ytickformat('%.2f');
xlim([0 max(ng_list)]);
legend(leg,'Location','northeast');
% ylim([0 3.2])
set(gca,'FontSize',20);
set(gca,'XTick',[0 1500 3000]);
set(gca,'YTick',[minneglogposts,(minneglogposts+maxneglogposts)/2,maxneglogposts]);
ylabel('-log{\it L}({\it n}_G|\lambda,\sigma)','FontSize',25);
xlabel('{\it n}_G','FontSize',25);
% title({'Method','Comparison'},'FontSize',25);
title('Posterior Minimization','FontSize',25);

if savenclose
    print('Figure_2c_zeiselposterior','-dtiffn')
    close
end

end