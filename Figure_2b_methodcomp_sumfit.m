function Figure_2b_methodcomp_sumfit(savenclose, directory)

if nargin < 2
    directory = [cd filesep 'MatFiles'];
    if nargin < 1
        savenclose = 0;
    end
end

load([directory filesep 'methodcomp_output_final.mat'],'output');

colamd_idx = [1:29 31];

curgen_mrmrplus(:,1) = output.mRMR_plus.nGen; % mRMR_plus = MRx3
curgen_dbscan(:,1) = flipud(output.DBSCAN.nGen(colamd_idx));
curgen_entropy(:,1) = output.Entropy.nGen;
curgen_colamd(:,1) = output.COLAMD.nGen(colamd_idx);
curgen_mrmr(:,1) = output.mRMR.nGen(colamd_idx);
ngenmat = [curgen_dbscan curgen_entropy curgen_colamd curgen_mrmr curgen_mrmrplus];

f2 = figure; hold on;
set(f2,'Position',[0 0 500 425]); hold on;

methodnames = fieldnames(output);
for i = 1:length(methodnames)
    if strcmp(methodnames{i},'COLAMD') || strcmp(methodnames{i},'mRMR')
        cellopt_met(:,i) = output.(methodnames{i}).LinR.pv(colamd_idx,3)...
            + output.(methodnames{i}).LinR.sst(colamd_idx,3)...
            + output.(methodnames{i}).LinR.vip(colamd_idx,3)...
            + output.(methodnames{i}).tau(colamd_idx,1)...
            + output.(methodnames{i}).LinR.micro(colamd_idx,3);
    elseif strcmp(methodnames{i},'DBSCAN')
        cellopt_met(:,i) = flipud(output.(methodnames{i}).LinR.pv(colamd_idx,3))...
            + flipud(output.(methodnames{i}).LinR.sst(colamd_idx,3))...
            + flipud(output.(methodnames{i}).LinR.vip(colamd_idx,3))...
            + flipud(output.(methodnames{i}).tau(colamd_idx,1))...
            + flipud(output.(methodnames{i}).LinR.micro(colamd_idx,3));
    else
        cellopt_met(:,i) = output.(methodnames{i}).LinR.pv(:,3)...
            + output.(methodnames{i}).LinR.sst(:,3)...
            + output.(methodnames{i}).LinR.vip(:,3)...
            + output.(methodnames{i}).tau...
            + output.(methodnames{i}).LinR.micro(:,3);
    end

end

figure(f2);
plot(ngenmat,cellopt_met,'LineWidth',2.5); hold on;
plot([540 540],[0 max(max(cellopt_met))+0.5],'k-','LineWidth',2); hold on;
ylim([min(min(cellopt_met))-0.5 max(max(cellopt_met))+0.5])
xlim([0 3855])
ylim([0 3.2])
set(gca,'FontSize',20);
set(gca,'XTick',[0 1750 3500]);
set(gca,'YTick',[0 1.1 2.3 3.5]);
ylabel('\Sigma_{fit}','FontSize',25);
xlabel('{\it n}_G','FontSize',25);
title('\Sigma_{fit} Across Subsetting Methods','FontSize',25);
legend({'DBSCAN','Entropy','colAMD','mRMR','MRx3'},'Location','southeast');

if savenclose
    print('methodcomp_sumfit','-dtiffn')
    close
end
end