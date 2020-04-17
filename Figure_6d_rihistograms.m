function Figure_6d_rihistograms(randstruct,savenclose)

if nargin < 2
    savenclose = 0;
end

for n = 1:length(randstruct.nclust)
    figure('Position',[0 0 650 550]); hold on;
    hist(randstruct.RandIndexRandom(n,:),20)
    plot([randstruct.RandIndexT(n) randstruct.RandIndexT(n)], [0 2200], 'r-', 'LineWidth', 3);
    ylim([0 2200]);
    xlabel('Rand Index','FontSize',60); 
    ylabel('No. Simulations','FontSize',60);
    title(sprintf('ARI = %.2f',randstruct.AdjustedRand(n)),'FontSize',60);
    mri = min(randstruct.RandIndexRandom(n,:)) - 0.05;
    set(gca, 'FontSize', 40);
    set(gca,'YTick',[]);
    if min(randstruct.RandIndexRandom(n,:)) < randstruct.RandIndexT(n)
        xlim([mri randstruct.RandIndexT(n)+0.05]);
        set(gca,'XTick',linspace(mri,randstruct.RandIndexT(n)+0.05,3));
    else
        xlim([randT-0.05,max(randstruct.RandIndexRandom(n,:))+0.05])
        set(gca,'XTick',linspace(randstruct.RandIndexT(n)-0.05,max(randstruct.RandIndexRandom(n,:))+0.05,3));
    end
    xtickformat('%.2f');
    box on;
    
    if savenclose
        print(sprintf('randindex_histogram_k_%d',randstruct.nclust(n)),'-dtiffn')
        close
    end
end
   
end
