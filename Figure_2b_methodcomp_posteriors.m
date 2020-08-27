function Figure_2b_methodcomp_posteriors(sigma, savenclose, directory)

if nargin < 3
    directory = [cd filesep 'MatFiles']; 
    if nargin < 2
        savenclose = 0;
        if nargin < 1
            sigma = 4400;
        end
    end
end
methodnames = {'colAMD', 'DBSCAN', 'Entropy', 'MRx3'};
loadstruct = load([directory filesep 'lambda250_gmm_cv5.mat'],'outstruct');
ng_mrx3 = zeros(1,length(loadstruct.outstruct));
neglogposts_mrx3 = ng_mrx3;
for i = 1:length(loadstruct.outstruct)
    ng_mrx3(i) = loadstruct.outstruct(i).nGen;
    negloglike = -(loadstruct.outstruct(i).likelihood^2); % outstruct.likelihood = GMM posterior
    neglogprior = (ng_mrx3(i)^2)/(2*sigma^2);
    neglogposts_mrx3(i) = negloglike + neglogprior;
end

loadstruct = load([directory filesep 'tasic_colAMD.mat'],'fitstruct');
ng_colamd = zeros(1,length(loadstruct.fitstruct));
neglogposts_colamd = ng_colamd;
for i = 1:length(loadstruct.fitstruct)
    ng_colamd(i) = loadstruct.fitstruct(i).nGen;
    negloglike = loadstruct.fitstruct(i).negloglikelihood;
    neglogprior = (ng_colamd(i)^2)/(2*sigma^2);
    neglogposts_colamd(i) = negloglike + neglogprior;
end

loadstruct = load([directory filesep 'tasic_DBSCAN.mat'],'fitstruct');
ng_dbscan = zeros(1,length(loadstruct.fitstruct));
neglogposts_dbscan = ng_dbscan;
for i = 1:length(loadstruct.fitstruct)
    ng_dbscan(i) = loadstruct.fitstruct(i).nGen;
    negloglike = loadstruct.fitstruct(i).negloglikelihood;
    neglogprior = (ng_dbscan(i)^2)/(2*sigma^2);
    neglogposts_dbscan(i) = negloglike + neglogprior;
end

loadstruct = load([directory filesep 'tasic_Entropy.mat'],'fitstruct');
ng_ent = zeros(1,length(loadstruct.fitstruct));
neglogposts_ent = ng_ent;
for i = 1:length(loadstruct.fitstruct)
    ng_ent(i) = loadstruct.fitstruct(i).nGen;
    negloglike = loadstruct.fitstruct(i).negloglikelihood;
    neglogprior = (ng_ent(i)^2)/(2*sigma^2);
    neglogposts_ent(i) = negloglike + neglogprior;
end

ngcell = cell(1,4); neglogpostcell = ngcell;
ngcell{4} = ng_mrx3; ngcell{1} = ng_colamd;
ngcell{2} = flipud(ng_dbscan); ngcell{3} = ng_ent;
neglogpostcell{4} = neglogposts_mrx3; neglogpostcell{1} = neglogposts_colamd;
neglogpostcell{2} = neglogposts_dbscan; neglogpostcell{3} = neglogposts_ent;

minneglogposts = zeros(1,length(neglogpostcell));
maxneglogposts = minneglogposts;
for i = 1:length(neglogpostcell)
    minneglogposts(i) = min(neglogpostcell{i});
    maxneglogposts(i) = max(neglogpostcell{i});
end

[~,minind] = min(minneglogposts);
bestmethod = methodnames{minind};
bestng = ngcell{minind}(neglogpostcell{minind} == minneglogposts(minind));

f2 = figure; hold on;
set(f2,'Position',[0 0 500 425]); hold on;
figure(f2);
for i = 1:length(methodnames)
    plot(ngcell{i},neglogpostcell{i},'LineWidth',2.5); hold on;
end
plot([540 540],[1.1*min(minneglogposts), 0.9*max(maxneglogposts)],'k-','LineWidth',2); hold on;
ylim([1.1*min(minneglogposts), 0.9*max(maxneglogposts)]); ytickformat('%.2f');
xlim([0 3855]);
% ylim([0 3.2])
set(gca,'FontSize',20);
set(gca,'XTick',[0 1750 3500]);
set(gca,'YTick',[min(minneglogposts), (min(minneglogposts)+max(maxneglogposts))/2,...
    max(maxneglogposts)]);
ylabel('-log{\it p}({\it n}_G|C)','FontSize',25);
xlabel('{\it n}_G','FontSize',25);
title('Posterior Probabilities','FontSize',25);
legend({'colAMD','DBSCAN','Entropy','MRx3',...
    sprintf(['Best nG = %d (' bestmethod ')'],bestng)},'Location','northeast');

if savenclose
    print('methodcomp_sumfit','-dtiffn')
    close
end
end