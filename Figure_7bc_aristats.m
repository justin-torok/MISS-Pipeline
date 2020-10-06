function Figure_7bc_aristats(randstruct,savenclose)

if nargin < 2
    savenclose = 0;
end

figure;
plot(1:length(randstruct.nclust),randstruct.StdAboveMean,'b-','LineWidth',3.5);
ylabel('SD Above Mean','FontSize',25);
if all(randstruct.StdAboveMean > 5)
    ylim([5,max(randstruct.StdAboveMean)+1])
    set(gca,'YTick',[6 13 20]);
end
xlabel('Number of Clusters','FontSize',25);
xlim([0.5 5.5]);
for i = 1:length(randstruct.nclust)
    xticklabs{i} = ['k = ' num2str(randstruct.nclust(i))];
end
set(gca,'XTickLabels',xticklabs);
set(gca,'FontSize',20);
title('MISS Rand Index vs. Null', 'FontSize', 30);

if savenclose
    print('ari_stdlineplot','-dtiffn');
    close
end

figure; 
bar_color = [1 0 0];
bar(1:length(randstruct.nclust),randstruct.AdjustedRand,0.75,'FaceColor',bar_color,'EdgeColor','k','LineWidth',1);
if all(randstruct.AdjustedRand > 0)
    ylim([0,1]);
    set(gca,'YTick',[0 0.5 1]);
end
ylabel('ARI','FontSize',25);
xlabel('Number of Clusters','FontSize',25);
xlim([0.5 5.5]);
for i = 1:length(randstruct.nclust)
    xticklabs{i} = ['k = ' num2str(randstruct.nclust(i))];
end
set(gca,'XTickLabels',xticklabs);
set(gca,'FontSize',20);
title('ARI Across Cluster Number', 'FontSize', 30);
if savenclose
    print('Figure_7_ari_barplot','-dtiffn');
    close
end

end
