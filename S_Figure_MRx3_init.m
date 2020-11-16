function S_Figure_MRx3_init(lambda,savenclose,directory)
if nargin < 3
    directory = [cd filesep 'MatFiles']; 
    if nargin < 2
        savenclose = 0;
        if nargin < 1
            lambda = 250;
        end
    end
end

% directory = '/Users/justintorok/Documents/MATLAB/CellType/MatFiles/MISS';
% lambda = 250;
load([directory filesep 'Tasic_Inputs.mat'],'genevct','voxvgene');
seed_true = 30;
geneinds_true = MRx3_Selector(genevct,voxvgene,3855,lambda,seed_true);

seeds_test = [10,20,50,80];
geneinds_test = zeros(length(seeds_test),3855);
for i = 1:length(seeds_test)
    geneinds_test(i,:) = MRx3_Selector(genevct,voxvgene,3855,lambda,seeds_test(i));
end

testngs = 50:3855;
overlapmat = zeros(length(seeds_test),length(testngs));
for i = 1:length(seeds_test)
    for j = 1:length(testngs)
        gis_true = geneinds_true(1:testngs(j));
        gis_test = geneinds_test(i,1:testngs(j));
        overlapmat(i,j) = sum(ismember(gis_true,gis_test))/testngs(j);
    end
end
lcell = cell(1,length(seeds_test));
figure; hold on;
for i = 1:length(seeds_test)
    plot(testngs,overlapmat(i,:),'LineWidth',2);
    lcell{i} = sprintf('mRMR Seed Size = %d',seeds_test(i));
end
xlabel('{\it n}_G'); ylabel('Overlap'); legend(lcell);
title({'Proportion Overlap Between MRx3' 'with Different Initialization Seed Sizes'});
set(gca,'FontSize',20);

if savenclose
    print('sfigmrx3init','-dtiffn')
end

end