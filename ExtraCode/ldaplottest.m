load('/Users/justintorok/Documents/MATLAB/CellType/MatFiles/MISS/lambda150_gmm_cv5.mat','outstruct');
for i = 1:length(outstruct)
    testnG(i) = outstruct(i).nGen;
    error_vec(i) = mean(outstruct(i).error);
    invposterior_vec(i) = 1 - mean(outstruct(i).posteriors);
end


nG_max = 3855;
nG_rescaled  = (testnG - min(testnG))/(nG_max - min(testnG));
error_rescaled = (error_vec - min(error_vec))/max(error_vec - min(error_vec));
sqdist = nG_rescaled.^2 + error_rescaled.^2;
[~,optind] = min(sqdist);
% optind = findchangepts(error_vec);
% optind2 = findchangepts(error_vec,'Statistic','std');
invpost_rescaled = (invposterior_vec - min(invposterior_vec))/max(invposterior_vec - min(invposterior_vec));
% optind3 = findchangepts(invposterior_vec);
% optind4 = findchangepts(invposterior_vec,'Statistic','std');
sqdist2 = nG_rescaled.^2 + invpost_rescaled.^2;
[~,optind2] = min(sqdist2);

figure('Units','Inches','Position',[0 0 20 10]);
subplot(1,2,1)
plot([testnG(optind), testnG(optind)],[0,max(error_vec)],'b','LineStyle','--');
hold on;
% plot([testnG(optind2), testnG(optind2)],[0,max(error_vec)],'g','LineStyle','--');
plot(testnG,error_vec,'Color',[1 0 0],'LineWidth',3);
hold off;
xlabel('nG'); ylabel('GMM Error');
legend({sprintf('Elbow Optimal - nG = %d', testnG(optind))});
title('GMM Error (5-fold CV) vs. nG');
set(gca,'FontSize',20);
subplot(1,2,2)
% plot([testnG(optind), testnG(optind)],[0,max(invposterior_vec)],'b','LineStyle','--');
hold on;
plot([testnG(optind2), testnG(optind2)],[0,max(invposterior_vec)],'g','LineStyle','--');
plot(testnG,invposterior_vec,'Color',[1 0 0],'LineWidth',3);
hold off;
xlabel('nG'); ylabel('1 - Posterior');
legend({sprintf('Elbow Optimal - nG = %d', testnG(optind2))});
title('GMM Posterior (5-fold CV) vs. nG');
set(gca,'FontSize',20);