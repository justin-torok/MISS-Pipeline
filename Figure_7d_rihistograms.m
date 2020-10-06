function Figure_7d_rihistograms(randstruct,savenclose,directory)

if nargin < 3
    directory = [cd filesep 'MatFiles'];
    if nargin < 2
        savenclose = 0;
        if nargin < 1
            randstruct = [];
        end
    end
end

if ~isempty(randstruct)
    for n = 1:length(randstruct.nclust)
        figure('Position',[0 0 650 550]); hold on;
        histogram(randstruct.RandIndexRandom(n,:),20)
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
else
    load([directory filesep 'tasic_l250_ng529.mat'],'outstruct');
    randstruct = Rand_Index_Calc(outstruct,1,'fore',directory);
    tasic_randT = randstruct.RandIndexT;
    randomT = randstruct.RandIndexRandom;
    load([directory filesep 'zeisel_l250_ng1168.mat'],'outstruct');
    randstruct = Rand_Index_Calc(outstruct,1,'fore',directory);
    zeisel_randT = randstruct.RandIndexT;
    for n = 1:length(randstruct.nclust)
        figure('Position',[0 0 650 550]); hold on;
        h = histogram(randomT(n,:),20); h.FaceColor = [0 1 0]; h.EdgeColor = [0 0 0];
        plot([tasic_randT(n) tasic_randT(n)], [0 2200], 'r-', 'LineWidth', 3);
        plot([zeisel_randT(n) zeisel_randT(n)], [0 2200], 'b-', 'LineWidth', 3);
        ylim([0 2200]);
        xlabel('Rand Index','FontSize',60); 
        ylabel('No. Simulations','FontSize',60);
        title(sprintf('k = %d',randstruct.nclust(n)),'FontSize',60);
        mri = min(randomT(n,:)) - 0.05;
        set(gca, 'FontSize', 40);
        set(gca,'YTick',[]);
        if min(randomT(n,:)) < min([zeisel_randT(n),tasic_randT(n)])
            xlim([mri max([zeisel_randT(n),tasic_randT(n)])+0.05]);
            set(gca,'XTick',linspace(mri,max([zeisel_randT(n),tasic_randT(n)])+0.05,3));
        else
            xlim([min([zeisel_randT(n),tasic_randT(n)])-0.05,max(randomT(n,:))+0.05])
            set(gca,'XTick',linspace(min([zeisel_randT(n),tasic_randT(n)])-0.05,max(randomT(n,:))+0.05,3));
        end
        xtickformat('%.2f');
        box on;

        if savenclose
            print(sprintf('randindex_histogram_k_%d',randstruct.nclust(n)),'-dtiffn')
            close
        end
    end    
end
end