function S_Figure_Residuals(nG,lambda,savenclose,directory)

if nargin < 4
    directory = [cd filesep 'MatFiles'];
    if nargin < 3
        savenclose = 0;
        if nargin < 2
            lambda = 250;
            if nargin < 1
                nG = 529;
            end
        end
    end
end
% lambda = 250; nG = 529; directory = '/Users/justintorok/Documents/MATLAB/CellType/MatFiles/MISS'; savenclose = 0;
load([directory filesep 'Tasic_Inputs.mat'],'genevct','voxvgene','gene_names','listBmap','nonzerovox');
sampleinds = [70,92,186];
listBvec = listBmap(nonzerovox);
[E_red, C_red] = GeneSelector(genevct,voxvgene,gene_names,nG,lambda,'MRx3');
D = CellDensityInference(E_red,C_red);
E_infer = C_red*D.';

residuals = zeros(1,size(E_infer,2));
for i = 1:length(residuals)
%     residuals(i) = sum(abs(E_infer(:,i) - E_red(:,i)))/sum(E_red(:,i));
    residuals(i) = sum(abs(E_infer(:,i) - E_red(:,i)));
end

samplebool = ismember(listBvec,sampleinds);
sample_errs = residuals(samplebool);
nonsample_errs = residuals(~samplebool);

datacell = {sample_errs,nonsample_errs};
rng('default');
xpos = zeros(length(residuals),1);
g = xpos;
xposscatter = @(x,y) 0.025 * (2*rand(length(x),1) - 1) + y;
plotindcell = {find(samplebool), find(~samplebool)};
for j = 1:length(plotindcell)
    curinds = plotindcell{j};
    g(curinds) = j;
    xpos(curinds) = xposscatter(curinds,j);
end

cmap = [[0.3 0.3 0.7];[0.7 0.3 0.3]];
f3 = figure('Position',[0 0 800 800]); hold on;
% v = gscatter(xpos,residuals,g,cmap,[],10,'off');
b = boxplot(residuals,g,'Colors',cmap,'Symbol','');
% v = violin(datacell,'facecolor',cmap);
set(b,{'linew'},{3})
set(gca, 'XTick', 1:2, 'XTickLabel', {'Sampled Regions','Unsampled Regions'});
% set(gca, 'YTick', [0,0.25,0.5,0.75,1])
% set(gca, 'XTickLabelRotation', 45);
ylabel('\Sigma_g |{\it E}_{red} - {\it C}_{red} * {\it D}|');
% ylabel('Rel. Error');
xlabel('');
set(gca,'TickLength',[0 0])
title('Sum of Absolute Error per Voxel');
set(gca, 'FontSize', 24);
if savenclose
    print('residuals','-dtiffn');
    close;
end

% % load pervoxLayerMaps.mat
% % load default_mousereg.mat
% 
load([directory filesep 'input_struct_voxelrender.mat'],'input_struct');
input_struct.voxUreg = 1;
input_struct.nbin = 5;
input_struct.centered = [1 2];
input_struct.xfac = 15;
input_struct.sphere = 1;
% input_struct.size_constant = 3;
input_struct.pointsize = 0.1;
input_struct.savenclose = savenclose;

% for j = 1:length(structIndex)
%     curres = residuals(listBvec==j);
%     residuals_reg(j) = mean(curres);
% end
% residuals_reg = ([residuals_reg(1:11) 0 residuals_reg(12:212) residuals_reg(1:11) 0 residuals_reg(12:212)]).';
testvals = residuals_reg;
% datamap = zeros(size(listBmap));
% for i = 1:length(structIndex)
%     [~,voxinds] = ismember(structIndex{i},nonzerovox);
%     allvox = nonzerovox(voxinds);
%     curvox = testvals(i);
%     datamap(allvox) = curvox;
%     clear allvox curvox voxinds
% end
% interpmap = datamap;
% interpmap = imresize3(interpmap,[133 81 115]);
% interpmap(interpmap<0) = 0;
% datinput = interpmap;
input_struct.data = testvals;
input_struct.cmap = [[0.3 0.3 0.7];[0.7 0.3 0.3]];
reggroups = 2*ones(426,1); reggroups([sampleinds+1,sampleinds+214]) = 1;
input_struct.region_groups = reggroups;
input_struct.img_labels = 'mean_residual_error';
figure;
brainframe(input_struct);
end

