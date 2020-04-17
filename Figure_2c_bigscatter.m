function Figure_2c_bigscatter(outstruct,idx,savenclose,directory)
% Designed to handle an outstruct of dimension length(ngenelist) that has a
% single field "bvals." idx is the index of interest that is within the
% range of ngenelist
if nargin < 4
    directory = [cd filesep 'MatFiles'];
    if nargin < 3
        savenclose = 0;
    end
end

ngenes = outstruct(idx).nGen;
neoinds = 57:94; 
foreinds = [1:11,23:56,141:156,170:212];
wbinds = [12:22 95:140 157:169];
indcell = {wbinds,foreinds,neoinds};
indcolor = {'b','m','g'};
indmark = {'*','o','d','s','^','p'};

LinRcalc = @(x,y) 2*corr(x,y)*std(x)*std(y)/(std(x)^2 + std(y)^2 + (mean(x) - mean(y))^2);

bigsim = 0; bigdat = 0;
f1 = figure; hold on;
for i = 1:length(indcell)
    load([directory filesep 'kim_totals_reorder_m.mat'],'kim_totals_reorder');
    testinds = indcell{i};
    kim_totals_reorder = kim_totals_reorder(testinds,:);
    nonnaninds = zeros(size(kim_totals_reorder));
    for j = 1:size(kim_totals_reorder,2)
        nonnaninds(:,j) = ~isnan(kim_totals_reorder(:,j));
    end
    nonnaninds = logical(nonnaninds);
    kim_pv_wb = kim_totals_reorder(nonnaninds(:,1),1);
    kim_sst_wb = kim_totals_reorder(nonnaninds(:,2),2);
    kim_vip_wb = kim_totals_reorder(nonnaninds(:,3),3);

    load([directory filesep 'murakami_totals_reorder.mat'],'murakami_totals_reorder');
    murakami_totals_reorder = murakami_totals_reorder(testinds);
    nonzeroinds = find(murakami_totals_reorder);
    murakami_tot = murakami_totals_reorder(nonzeroinds);

    sumcts_wb = outstruct(idx).Bsums(1:213,:) + outstruct(idx).Bsums(214:end,:);
    sumcts_wb(12,:) = [];
    sumcts_wb_mod = sumcts_wb(testinds,:);
    sumcts_pv_wb = sumcts_wb_mod(nonnaninds(:,1),3);
    sumcts_sst_wb = sumcts_wb_mod(nonnaninds(:,2),6);
    sumcts_vip_wb = sumcts_wb_mod(nonnaninds(:,3),7);
    sumcts_sum = sum(sumcts_wb_mod(nonzeroinds,:),2);
    
    load([directory filesep 'keller_micro_totals_reorder.mat'],'keller_micro_totals_reorder');
    load([directory filesep 'keller_micro_listB.mat'],'keller_micro_listB');
    load([directory filesep 'herculano_houzel_neuron_totals_reorder.mat'],'hh_neuron_totals_reorder');
    load([directory filesep 'herculano_houzel_neuron_listB.mat'],'herculano_houzel_neuron_listB');
    
    sumcts_micro_wb = sumcts_wb(:,23);
    sumcts_neuron_wb = sum(sumcts_wb(:,1:21),2);
    
    hh_micro = [];
    infer_micro = [];
    for j = 1:size(keller_micro_listB)
        hh = keller_micro_listB(j,:);
        hhinds = hh(hh > 0);
        if all(ismember(hhinds,testinds))
            hh_micro = [hh_micro, sum(keller_micro_totals_reorder(hhinds))];
            infer_micro = [infer_micro, sum(sumcts_micro_wb(hhinds))];
        end
    end
    
    hh_neuron = [];
    infer_neuron = [];
    for j = 1:size(herculano_houzel_neuron_listB)
        hh = herculano_houzel_neuron_listB(j,:);
        hhinds = hh(hh > 0);
        if all(ismember(hhinds,testinds))
            hh_neuron = [hh_neuron, sum(hh_neuron_totals_reorder(hhinds))];
            infer_neuron = [infer_neuron, sum(sumcts_neuron_wb(hhinds))];
        end
    end
    
    cellsim = [sumcts_pv_wb; sumcts_sst_wb; sumcts_vip_wb;...
        sumcts_sum; infer_micro.'; infer_neuron.'];
    
    celldat = [kim_pv_wb; kim_sst_wb; kim_vip_wb;...
        murakami_tot; hh_micro.'; hh_neuron.'];
    
    cellsim_cell = {sumcts_pv_wb, sumcts_sst_wb, sumcts_vip_wb,...
                    sumcts_sum, infer_micro.', infer_neuron.'};
                
    celldat_cell = {kim_pv_wb, kim_sst_wb, kim_vip_wb,...
                    murakami_tot, hh_micro.', hh_neuron.'};
    
    bigdat = [bigdat; celldat];
    bigsim = [bigsim; cellsim];
    
    figure(f1); hold on;
    for j = 1:length(indmark)
        if ismember(j,[1])
            scatter(cellsim_cell{j},celldat_cell{j},50,[indcolor{i} indmark{j}]);
        else
            scatter(cellsim_cell{j},celldat_cell{j},50,[indcolor{i} indmark{j}],'filled');
        end
    end

end

bigdat(1) = []; bigsim(1) = [];
figure(f1); hold on;
plot([1 max(cat(1,bigsim,bigdat))],[1 max(cat(1,bigsim,bigdat))],'k-'); hold on;
ylim([1 max(bigdat)]); xlim([1 max(bigdat)]);
yticks([10^0 10^2 10^4 10^6]);
xticks([10^0 10^2 10^4 10^6]);
set(gca,'XScale','log');
set(gca,'Yscale','log');
set(gca,'FontSize',18);
ylabel('Literature Counts','FontSize',20);
xlabel('Inferred Counts','FontSize',20);
title('MISS + MRx3, Literature vs. Inferred','FontSize',20);
txt = {['{\it n}_G = ' num2str(ngenes)]; sprintf('R_C = %.2f',LinRcalc(bigsim,bigdat));...
    sprintf('R = %.2f',corr(bigsim,bigdat))};
text(10,2.5*10^5,txt,'FontSize',15);

if savenclose
    print('alldata_scatter','-dtiffn');
    close
end
end