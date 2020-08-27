sigmas = [4350];
matdir = '/Users/justintorok/Documents/MATLAB/CellType/MatFiles';

load([matdir filesep sprintf('zeisel_gmm_l250_%d.mat',1)],'fitstruct');
for j = 2:4
    x = load([matdir filesep sprintf('zeisel_gmm_l250_%d.mat',j)],'fitstruct');
    fitstruct = cat(2,fitstruct,x.fitstruct);
end
    
for k = 1:length(sigmas)
    for i = 1:length(fitstruct)
        ng_param_list(i) = fitstruct(i).nGen;
        negloglikelihoods(i,k) = fitstruct(i).negloglikelihood;
        neglogpriors(i,k) = (fitstruct(i).nGen^2)/(2*sigmas(k)^2);        
    end
%     templeg{k} = sprintf('Sigma = %d', sigmas(k));
end

neglogposteriors = negloglikelihoods + neglogpriors;
mininds = zeros(1,length(sigmas));
for k = 1:size(mininds,2)
    [~,minind] = min(neglogposteriors(:,k));
    mininds(k) = minind;
end

figure('Units','Inches','Position',[0 0 20 10]);
cmap = hsv(length(sigmas));
for k = 1:length(sigmas)
%     plot(ng_param_list,squeeze(negloglikelihoods(:,j,1)),'Color',[0 0 0],'LineWidth',1);
    hold on;
    plot(ng_param_list,squeeze(neglogposteriors(:,k)),'Color',cmap(k,:),'LineWidth',2);
end
for k = 1:length(sigmas)
    scatter(ng_param_list(mininds(k)),neglogposteriors(mininds(k),k),250,cmap(k,:),'LineWidth',1.5);
%     templeg = cell(1,4);
%     templeg{1} = "-log(Likelihood)";
    templeg{k} = sprintf("-log(Posterior), sigma = %d - Minimum @ nG = %d",sigmas(k),ng_param_list(mininds(k)));
end
legend(templeg); xlabel('nG'); ylabel('-log(P)');
title(sprintf('GMM Posteriors, lambda = %d, no CV', fitstruct(1).lambda));
ylim([min(min(negloglikelihoods)),max(max(negloglikelihoods))])
set(gca,'FontSize',16);
    


