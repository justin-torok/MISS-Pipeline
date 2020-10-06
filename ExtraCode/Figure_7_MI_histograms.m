function Figure_7_MI_histograms(mistruct,savenclose)

if nargin < 2
    savenclose = 0;
end

for n = 1:length(mistruct.nclust)
    figure('Position',[0 0 650 550]); hold on;
    hist(mistruct.MutualInfoRandom(n,:),20)
    plot([mistruct.MutualInfoT(n) mistruct.MutualInfoT(n)], [0 2200], 'r-', 'LineWidth', 3);
    ylim([0 2200]);
    xlabel('Mutual Information','FontSize',60); 
    ylabel('No. Simulations','FontSize',60);
    title(sprintf('MI %sbrain, k = %d',mistruct.onttype,mistruct.nclust(n)),'FontSize',60);
    mmi = min(mistruct.MutualInfoRandom(n,:)) - 0.05;
    set(gca, 'FontSize', 40);
    set(gca,'YTick',[]);
    if min(mistruct.MutualInfoRandom(n,:)) < mistruct.MutualInfoT(n)
        xlim([mmi mistruct.MutualInfoT(n)+0.05]);
        set(gca,'XTick',linspace(mmi,mistruct.MutualInfoT(n)+0.05,3));
    else
        xlim([mistruct.MutualInfoT(n)-0.05,max(mistruct.MutualInfoRandom(n,:))+0.05])
        set(gca,'XTick',linspace(mistruct.MutualInfoT(n)-0.05,max(mistruct.MutualInfoRandom(n,:))+0.05,3));
    end
    xtickformat('%.2f');
    box on;
    
    if savenclose
        print(sprintf('MutualInfo_histogram_k_%d',mistruct.nclust(n)),'-dpng')
        close
    end
end
   
end
