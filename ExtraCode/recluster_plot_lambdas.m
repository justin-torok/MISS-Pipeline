matdir = '/Users/justintorok/Documents/MATLAB/CellType/MatFiles/MISS';
mats2load = {'lambda150_gmm_cv5.mat','lambda150_gmm_nocv.mat'};
elbow_crit = 'plateau'; %options: plateau+2std, prox2origin, plateau, below2.5pct, kneepoint, changepoint

for j = 1:length(mats2load)
    load([matdir filesep mats2load{j}],'outstruct','ng_param_list');
    for i = 1:length(outstruct)
        error(i,j) = 1-mean(outstruct(i).posteriors);
        accuracy(i,j) = outstruct(i).accuracy;
%         std(i,j) = outstruct(i).std;
    end
    lambda(j) = outstruct(1).lambda;
    clear outstruct
end

normnG = ((ng_param_list - min(ng_param_list)) / (3855 - min(ng_param_list))).';
colspecs = lines(length(mats2load));
figure; hold on;
for k = 1:length(mats2load)
    norm_error = (error(:,k) - min(error(:,k))) / (max(error(:,k)) - min(error(:,k)));
    if strcmp(elbow_crit,'prox2origin')
        dist2origin = sqrt((normnG-0).^2 + (norm_error-0).^2);
        [~,elbowind] = min(dist2origin);
        elbows(k) = ng_param_list(elbowind);
    elseif strcmp(elbow_crit,'plateau')
        end_error = error(end,k);
        at_below_end = find(error(:,k)<=end_error);
        elbowind = at_below_end(1);
        elbows(k) = ng_param_list(elbowind);
    elseif strcmp(elbow_crit,'plateau+2std')
        stdend2 = error(end,k) + 2 * std(end,k);
        at_below_stdend2 = find(error(:,k)<=stdend2);
        elbowind = at_below_stdend2(1);
        elbows(k) = ng_param_list(elbowind);
    elseif strcmp(elbow_crit,'below2.5pct')
        at_below_thresh = find(error(:,k)<=0.025);
        elbowind = at_below_thresh(1);
        elbows(k) = ng_param_list(elbowind);
    elseif strcmp(elbow_crit,'kneepoint')
        [~,elbowind] = knee_pt(norm_error);
        elbows(k) = ng_param_list(elbowind);
    elseif strcmp(elbow_crit,'changepoint')
        elbowind = findchangepts(error(:,k),'Statistic','mean');
        elbows(k) = ng_param_list(elbowind);
    end
    elbow_error(k) = error(elbowind,k);
    lambdalabs{2*k-1} = num2str(lambda(k));
    lambdalabs{2*k} = 'Elbow';
    plot(normnG,norm_error,'Color',colspecs(k,:),'LineWidth',2); hold on;
    plot([normnG(elbowind) normnG(elbowind)],[0 1],'Color',colspecs(k,:),'LineWidth',1,'LineStyle','--'); hold on;
end
xlabel('Norm. nG','FontSize',18);
ylabel('Norm. Error','FontSize',18);
xticks([0 0.2 0.4 0.6 0.8 1]);
yticks([0 0.25 0.5 0.75 1]);
set(gca,'FontSize',16);
title('NCC Elbow Curve','FontSize',18);
legend(lambdalabs,'Location','northeast');

figure;
plot(repmat(ng_param_list.',1,length(mats2load)),error,'LineWidth',2.5); hold on;
plot([elbows; elbows],repmat([0; max(max(error))],1,length(mats2load)),'LineWidth',1,'LineStyle','--');
yticks([0 0.02 0.04 0.06]);
ylabel('NCC Error','FontSize',18);
xticks([500 1000 1500 2000 2500 3000 3500]);
xlabel('nG','FontSize',18);
set(gca,'FontSize',16);
title('NCC Error vs. nG as a Function of \lambda','FontSize',18);
legend(lambdalabs,'Location','northeast');

mean_error = mean(error);
error_quant = [elbow_error; mean_error];
error_quant = error_quant.';
figure; hold on;
bar(error_quant,'LineWidth',0.000001); hold on;
xticks(1:7);
xticklabels({lambdalabs{1},lambdalabs{2}});
% xticklabels({lambdalabs{1},lambdalabs{3},lambdalabs{5},lambdalabs{7},...
%     lambdalabs{9},lambdalabs{11},lambdalabs{13}});
yticks([0.02 0.03]);
ylim([0.02 0.03]);
xlabel('\lambda','FontSize',18);
ylabel('NCC Error','FontSize',18);
set(gca,'FontSize',14);
title('Overall NCC Error Assessments Per \lambda','FontSize',18);
legend({'Error at Elbow','Mean Error'});


