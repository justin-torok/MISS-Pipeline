matdir = '/Users/justintorok/Documents/MATLAB/CellType/MatFiles/MISS';
mats2load = {'lambda0_gmm_cv5.mat','lambda50_gmm_cv5.mat',...
    'lambda100_gmm_cv5.mat','lambda150_gmm_cv5.mat',...
    'lambda200_gmm_cv5.mat','lambda250_gmm_cv5.mat',...
    'lambda300_gmm_cv5.mat','lambda350_gmm_cv5.mat'};
sigmas = [4000,4100,4200,4300,4400,4500,4600,4700,4800];
% sigmas = [3855, 5*3855, 10*3855, 15*3855, 20*3855, 25*3855];
cvs = {'5-Fold CV'};

for k = 1:length(sigmas)
    for j = 1:length(mats2load)
        load([matdir filesep mats2load{j}],'outstruct','ng_param_list');
        for i = 1:length(outstruct)
            error(i,j) = outstruct(i).error;
            accuracy(i,j) = outstruct(i).accuracy;
            if outstruct(i).lambda == 150
                negloglikelihoods(i,j,k) = -(outstruct(i).posteriors^2);
            else
                negloglikelihoods(i,j,k) = -(outstruct(i).likelihood^2);
            end
            neglogpriors(i,j,k) = (ng_param_list(i)^2)/(2*sigmas(k)^2);
        end
        lambda(j) = outstruct(1).lambda;
        clear outstruct
    end
    templeg{k} = sprintf('Sigma = %d', sigmas(k));
end
neglogposteriors = negloglikelihoods + neglogpriors;
mininds = zeros(length(mats2load),length(sigmas));
for j = 1:size(mininds,1)
    for k = 1:size(mininds,2)
        [~,minind] = min(neglogposteriors(:,j,k));
        mininds(j,k) = minind;
    end
end
figure('Units','Inches','Position',[0 0 25 15]);
cmap = hsv(length(lambda));
for k = 1:length(sigmas)
    subplot(sqrt(length(sigmas)),sqrt(length(sigmas)),k);
%     plot(ng_param_list,squeeze(negloglikelihoods(:,j,1)),'Color',[0 0 0],'LineWidth',1);
    hold on;
    for j = 1:length(mats2load)
        plot(ng_param_list,squeeze(neglogposteriors(:,j,k)),'Color',cmap(j,:),'LineWidth',2);
    end
    for j = 1:length(mats2load)
        scatter(ng_param_list(mininds(j,k)),neglogposteriors(mininds(j,k),j,k),250,cmap(j,:),'LineWidth',1.5);
    end
    templeg = cell(1,length(mats2load));
%     templeg{1} = "-log(Likelihood)";
    for j = 1:length(mats2load)
        templeg{j} = sprintf("-log(Posterior), lambda = %d - Minimum @ nG = %d", lambda(j),ng_param_list(mininds(j,k)));
    end
    title(sprintf(['GMM Posteriors, sigma = %d, ' cvs{1}], sigmas(k)));
    xlabel('nG'); ylabel('-log(P)');
    ylim([min(min(negloglikelihoods(:,j,:))),max(max(negloglikelihoods(:,j,:)))])
    set(gca,'FontSize',16);
    legend(templeg,'FontSize',12); 
end


% 
% normnG = ((ng_param_list - min(ng_param_list)) / (3855 - min(ng_param_list))).';
% colspecs = lines(length(mats2load));
% figure; hold on;
% for k = 1:length(mats2load)
%     norm_error = (error(:,k) - min(error(:,k))) / (max(error(:,k)) - min(error(:,k)));
%     if strcmp(elbow_crit,'prox2origin')
%         dist2origin = sqrt((normnG-0).^2 + (norm_error-0).^2);
%         [~,elbowind] = min(dist2origin);
%         elbows(k) = ng_param_list(elbowind);
%     elseif strcmp(elbow_crit,'plateau')
%         end_error = error(end,k);
%         at_below_end = find(error(:,k)<=end_error);
%         elbowind = at_below_end(1);
%         elbows(k) = ng_param_list(elbowind);
%     elseif strcmp(elbow_crit,'plateau+2std')
%         stdend2 = error(end,k) + 2 * std(end,k);
%         at_below_stdend2 = find(error(:,k)<=stdend2);
%         elbowind = at_below_stdend2(1);
%         elbows(k) = ng_param_list(elbowind);
%     elseif strcmp(elbow_crit,'below2.5pct')
%         at_below_thresh = find(error(:,k)<=0.025);
%         elbowind = at_below_thresh(1);
%         elbows(k) = ng_param_list(elbowind);
%     elseif strcmp(elbow_crit,'kneepoint')
%         [~,elbowind] = knee_pt(norm_error);
%         elbows(k) = ng_param_list(elbowind);
%     elseif strcmp(elbow_crit,'changepoint')
%         elbowind = findchangepts(error(:,k),'Statistic','mean');
%         elbows(k) = ng_param_list(elbowind);
%     end
%     elbow_error(k) = error(elbowind,k);
%     lambdalabs{2*k-1} = num2str(lambda(k));
%     lambdalabs{2*k} = 'Elbow';
%     plot(normnG,norm_error,'Color',colspecs(k,:),'LineWidth',2); hold on;
%     plot([normnG(elbowind) normnG(elbowind)],[0 1],'Color',colspecs(k,:),'LineWidth',1,'LineStyle','--'); hold on;
% end
% xlabel('Norm. nG','FontSize',18);
% ylabel('Norm. Error','FontSize',18);
% xticks([0 0.2 0.4 0.6 0.8 1]);
% yticks([0 0.25 0.5 0.75 1]);
% set(gca,'FontSize',16);
% title('NCC Elbow Curve','FontSize',18);
% legend(lambdalabs,'Location','northeast');
% 
% figure;
% plot(repmat(ng_param_list.',1,length(mats2load)),error,'LineWidth',2.5); hold on;
% plot([elbows; elbows],repmat([0; max(max(error))],1,length(mats2load)),'LineWidth',1,'LineStyle','--');
% yticks([0 0.02 0.04 0.06]);
% ylabel('NCC Error','FontSize',18);
% xticks([500 1000 1500 2000 2500 3000 3500]);
% xlabel('nG','FontSize',18);
% set(gca,'FontSize',16);
% title('NCC Error vs. nG as a Function of \lambda','FontSize',18);
% legend(lambdalabs,'Location','northeast');
% 
% mean_error = mean(error);
% error_quant = [elbow_error; mean_error];
% error_quant = error_quant.';
% figure; hold on;
% bar(error_quant,'LineWidth',0.000001); hold on;
% xticks(1:7);
% xticklabels({lambdalabs{1},lambdalabs{2}});
% % xticklabels({lambdalabs{1},lambdalabs{3},lambdalabs{5},lambdalabs{7},...
% %     lambdalabs{9},lambdalabs{11},lambdalabs{13}});
% yticks([0.02 0.03]);
% ylim([0.02 0.03]);
% xlabel('\lambda','FontSize',18);
% ylabel('NCC Error','FontSize',18);
% set(gca,'FontSize',14);
% title('Overall NCC Error Assessments Per \lambda','FontSize',18);
% legend({'Error at Elbow','Mean Error'});
% 
% 
