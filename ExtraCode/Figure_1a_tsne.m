function Figure_1a_tsne(savenclose, directory)
% Note - result may look different from the Figure 1a panel in the
% manuscript because a default rng seed was not used to generate it.

if nargin < 2
    directory = [cd filesep 'MatFiles'];
    if nargin < 1
        savenclose = 0;
    end
end

load([directory filesep 'PresetInputs.mat'],'meanexprmat');
C = meanexprmat;
ctmean = mean(C,1);
ctmean = repmat(ctmean,size(C,1),1);
C = C ./ ctmean;
Y = tsne(C);
[idx] = kmeans(C,25);
figure;
gscatter(Y(:,1), Y(:,2),idx,[],[],10,0)
xlabel('tSNE Dimension 1'); ylabel('tSNE Dimension 2');
set(gca,'XTickLabel',[],'YTickLabel',[]);
set(gca,'FontSize',20);
if savenclose
    saveas(gcf,'tSNE_C','tiffn');
    close
end
end