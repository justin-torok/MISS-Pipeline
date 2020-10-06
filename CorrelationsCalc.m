function [LinRvals,PeaRvals,SpearRhovals] = CorrelationsCalc(outstruct,idx,directory)
% Calculates Lin and Pearson correlations from a pre-created outstruct,
% with a specified index (idx) that indicates which row of the outstruct to
% calculate these metrics for.

if nargin < 3
    directory = [cd filesep 'MatFiles'];
end

% Define gross anatomical regions of interest
amg = 1:11; sub = 23:25; hip = 26:36; hyp = 37:56; neo = 57:94; 
olf = 141:148; pal = 149:156; str = 170:177; tha = 178:212;
cer = 12:22; med = 95:119; mid = 120:140; pns = 157:169;

neoinds = neo; % neocortical regions in the CCF atlas
foreinds = [amg sub hip hyp olf pal str tha]; % forebrain regions in the CCF atlas minus neocortex
hindinds = [cer med mid pns]; % mid/hindbrain regions in the CCF atlas
wbinds = 1:212; % all 212 (bilateral) regions in the CCF atlas
indcell = {neoinds,foreinds,hindinds,wbinds};

% The following function calculates the Lin's Concordance Correlation
% Coefficient.
LinRcalc = @(x,y) 2*corr(x,y)*std(x)*std(y)/(std(x)^2 + std(y)^2 + (mean(x) - mean(y))^2);

for i = 1:length(indcell)
    testinds = indcell{i};
    
    % Extract Kim et al, 2017 interneuron data from pre-created file
    load([directory filesep 'kim_totals_reorder_m.mat']);
    kim_totals_reorder = kim_totals_reorder(testinds,:);
    nonnaninds = zeros(size(kim_totals_reorder));
    for j = 1:size(kim_totals_reorder,2)
        nonnaninds(:,j) = ~isnan(kim_totals_reorder(:,j));
    end
    nonnaninds = logical(nonnaninds);
    kim_pv_wb = kim_totals_reorder(nonnaninds(:,1),1);
    kim_sst_wb = kim_totals_reorder(nonnaninds(:,2),2);
    kim_vip_wb = kim_totals_reorder(nonnaninds(:,3),3);
    
    % Extract Murakami et al, 2018 total cell data from pre-created file
    load([directory filesep 'murakami_totals_reorder.mat']); 
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
    
    % Extract Keller et al, 2018 microglia data from pre-created files
    load([directory filesep 'keller_micro_totals_reorder.mat']);
    load([directory filesep 'keller_micro_listB.mat']);
    
    sumcts_micro_wb = sumcts_wb(:,23);
    sumcts_neuron_wb = sum(sumcts_wb(:,1:21),2);

    keller_micro = [];
    infer_micro = [];
    for j = 1:size(keller_micro_listB,1)
        kk = keller_micro_listB(j,:);
        kkinds = kk(kk > 0);
        if all(ismember(kkinds,testinds))
            keller_micro = [keller_micro, sum(keller_micro_totals_reorder(kkinds))];
            infer_micro = [infer_micro, sum(sumcts_micro_wb(kkinds))];
        end
    end
    
    % Extract Herculano-Houzel et al, 2013 total neuron data from
    % pre-created files
    load([directory filesep 'herculano_houzel_neuron_totals_reorder.mat']);
    load([directory filesep 'herculano_houzel_neuron_listB.mat']);
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
    
    LinRvals.pv(i) = LinRcalc(kim_pv_wb,sumcts_pv_wb);
    LinRvals.sst(i) = LinRcalc(kim_sst_wb,sumcts_sst_wb);
    LinRvals.vip(i) = LinRcalc(kim_vip_wb,sumcts_vip_wb);
    LinRvals.all(i) = LinRcalc(murakami_tot,sumcts_sum);
    LinRvals.micro(i) = LinRcalc(keller_micro.',infer_micro.');
%     LinRvals.neuron(i) = LinRcalc(hh_neuron.',infer_neuron.');
    
    PeaRvals.pv(i) = corr(kim_pv_wb,sumcts_pv_wb);
    PeaRvals.sst(i) = corr(kim_sst_wb,sumcts_sst_wb);
    PeaRvals.vip(i) = corr(kim_vip_wb,sumcts_vip_wb);
    PeaRvals.all(i) = corr(murakami_tot,sumcts_sum);
    PeaRvals.micro_pearson(i) = corr(keller_micro.',infer_micro.');
%     PeaRvals.neuron_pearson(i) = corr(hh_neuron.',infer_neuron.');
    
    SpearRhovals.pv(i) = corr(kim_pv_wb,sumcts_pv_wb,'type','Spearman');
    SpearRhovals.sst(i) = corr(kim_sst_wb,sumcts_sst_wb,'type','Spearman');
    SpearRhovals.vip(i) = corr(kim_vip_wb,sumcts_vip_wb,'type','Spearman');
    SpearRhovals.all(i) = corr(murakami_tot,sumcts_sum,'type','Spearman');
    SpearRhovals.micro_pearson(i) = corr(keller_micro.',infer_micro.','type','Spearman');
%     SpearRhovals.neuron_pearson(i) = corr(hh_neuron.',infer_neuron.','type','Spearman');
end