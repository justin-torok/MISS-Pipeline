% types = {'Pvalb','Sst','Vip','L4','L5_IT','Oligo','Endo'};
types = [9:15];
% load mrmrplus_lambda150_finerange.mat;
load mrmrplus_540_150_outstruct.mat
method = 'mRMRplus';
peakind = 1;

% load sordes_output_allrange.mat;
% method = 'SoRDES';
% peakind = 4;
load pervoxLayerMaps.mat
load default_mousereg.mat
input_struct.voxUreg = 0;
input_struct.nbin = 5;
input_struct.centered(1) = 0;
input_struct.xfac = 2;
input_struct.pointsize = 0.1;
input_struct.savenclose = 0;
% VoxelRender_3D_pccs(outstruct,peakind,method,types,1);
for j = 3
%     testvals = outstruct(peakind).Bvals(:,types(j));
    testvals = outstruct.Bvals(:,types(j));
    datamap = zeros(size(GENGDmod));
    for i = 1:length(structList)
        [~,voxinds] = ismember(structIndex{i},nonzerovox);
        allvox = nonzerovox(voxinds);
        curvox = testvals(voxinds);
        datamap(allvox) = curvox;
        clear allvox curvox voxinds
    end
    interpmap = datamap;
    interpmap = imresize3(interpmap,[133 81 115]);
    interpmap(interpmap<0) = 0;
    datinput = interpmap;
    colmap = [ones(5,1) flipud((0.2:0.2:1)') zeros(5,1)];
    input_struct.data = datinput;
    input_struct.cmap = colmap;
    input_struct.img_labels = classkey{2,types(j)};
%     figure;
%     brainframe(input_struct);
end

regsums = sum(outstruct.Bsums);
regsums = repmat(regsums,426,1);
ct_regprop = outstruct.Bsums ./ regsums;
ct_regprop([12 12+213],:) = [];
amg = 1:11; cer = 12:22; sub = 23:25; hip = 26:36; hyp = 37:56; neo = 57:94;
med = 95:119; mid = 120:140; olf = 141:148; pal = 149:156; pns = 157:169;
str = 170:177; tha = 178:212;
sumprop(1,:) = sum(ct_regprop(amg,:));
sumprop(2,:) = sum(ct_regprop(cer,:));
sumprop(3,:) = sum(ct_regprop(sub,:));
sumprop(4,:) = sum(ct_regprop(hip,:));
sumprop(5,:) = sum(ct_regprop(hyp,:));
sumprop(6,:) = sum(ct_regprop(neo,:));
sumprop(7,:) = sum(ct_regprop(med,:));
sumprop(8,:) = sum(ct_regprop(mid,:));
sumprop(9,:) = sum(ct_regprop(olf,:));
sumprop(10,:) = sum(ct_regprop(pal,:));
sumprop(11,:) = sum(ct_regprop(pns,:));
sumprop(12,:) = sum(ct_regprop(str,:));
sumprop(13,:) = sum(ct_regprop(tha,:));
labkey = classkey(2,:);
labkey{9} = 'L2+3-IT';
labeky{11} = 'L5-IT';
labkey{12} = 'L5-PT';
labkey{13} = 'L6-CT';
labkey{14} = 'L6-IT';
MajRegNms = {'Amg','Cer','Sub','Hip','Hyp','Neo','Med','Mid','Olf','Pal','Pns','Str','Tha'};

figure('Position',[0 0 600 150]);
imagesc(sumprop(end,:),[0.01 0.1]); colormap copper;
xticks(1:25);
xticklabels(labkey);
xtickangle(90);
yticks([]);
set(gca,'TickLength',[0 0]);
% yticks(1:13);
% yticklabels(MajRegNms);
set(gca,'FontSize',16);
title('CT Concentration in Thalamus','FontSize',18);


regtotals = outstruct.Bsums;
normtotals = zeros(426,25);
regentropy = zeros(1,426);
entropycalc = @(x) -dot(log(x+eps),x);
for i = 1:426
   normtotals(i,:) = regtotals(i,:)/sum(regtotals(i,:)); 
   regentropy(i) = entropycalc(normtotals(i,:));
end
regentropy(isnan(regentropy)) = 0;
datinput = ones(426,1)./regentropy.';
datinput([12 12+213]) = [];

newinput = (datinput(1:212) + datinput(213:end)) / 2;
reglabs = {'Amg','Cer','Sub','Hip','Hyp','Neo','Med','Mid','Olf','Pal','Pns','Str','Tha';...
            1:11,12:22,23:25,26:36,37:56,57:94,95:119,120:140,141:148,149:156,157:169,170:177,178:212};
plotlabs = reglabs(1,:);
rng('default');
xpos = zeros(212,1);
g = xpos;
xposscatter = @(x,y) 0.2 * (2*rand(length(x),1) - 1) + y;
for j = 1:length(plotlabs)
    curinds = reglabs{2,j};
    g(curinds) = j;
    xpos(curinds) = xposscatter(curinds,j);
%     summaryvec = zeros(2,1);
%     summaryvec(1) = mean(neuronratio);
%     summaryvec(2) = std(neuronratio);
%     summarymat(:,j) = summaryvec;
end

cmap = hsv(13);
f3 = figure('Position',[0 0 500 1000]); hold on; box on;
v = gscatter(newinput,xpos,g,cmap,[],20,'off');
b = boxplot(newinput,g,'orientation','horizontal','Colors',cmap,'Symbol','');
set(b,{'linew'},{3})
% for m = 1:size(plotlabs)
%     val1 = summarymat(1,m) + summarymat(2,m);
%     val2 = summarymat(1,m) - summarymat(2,m);
%     plot([m-0.3, m+0.3], [val1 val1], ':', 'Color', cmap(m,:), 'LineWidth', 2)
%     plot([m-0.3, m+0.3], [val2 val2], ':', 'Color', cmap(m,:), 'LineWidth', 2)
%     plot([m-0.3, m+0.3], [summarymat(1,m) summarymat(1,m)], '-', 'Color', cmap(m,:), 'LineWidth', 3)
% end
% plot([0.5,0.5],[0,14],'k--','LineWidth',2);
plot([median(newinput),median(newinput)],[0,14],'k--','LineWidth',2);
set(gca, 'YTick', 1:13, 'YTickLabel', reglabs(1,:));
set(gca, 'XTick', [0.5 2 3.5])
% set(gca, 'XTickLabelRotation', 45);
% ylabel('Proportion Glutamatergic Neurons');
ylabel('');
xlabel('');
xlim([0 max(datinput)]);
set(gca,'TickLength',[0 0])
title({sprintf('Homogeneity Index \n Across Regions')});
set(gca, 'FontSize', 30,'LineWidth',0.75);
% view([90 -90])


