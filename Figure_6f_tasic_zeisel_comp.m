function Figure_6f_tasic_zeisel_comp(savenclose,directory)

if nargin < 2
    directory = [cd filesep 'MatFiles'];
    if nargin < 1
        savenclose = 0;
    end
end

LinRcalc = @(x,y) 2*corr(x,y)*std(x)*std(y)/(std(x)^2 + std(y)^2 + (mean(x) - mean(y))^2);
tasic = load([directory filesep 'tasic_l250_ng529.mat'],'outstruct');
Bcorrected = tasic.outstruct.corrB;
[sumB,meanB] = Voxel_To_Region(Bcorrected,directory);
tasic.outstruct.Bsums = sumB;
tasic.outstruct.Bmeans = meanB;

zeisel = load([directory filesep 'zeisel_l250_nG1168_outputs.mat'],'outstruct');
Bcorrected = zeisel.outstruct.corrB;
[sumB,meanB] = Voxel_To_Region(Bcorrected,directory);
zeisel.outstruct.Bsums = sumB;
zeisel.outstruct.Bmeans = meanB;


ctypes.Zeisel = {[150:173 35:36 21:25 47:56 84:95 9 13 18];...%Glut
    [174:194 26:33 96:109 57:65 128:132 118:123 136:137 10:11 14];...%GABA
    [16 17 124 125 110 111 115 116 117];...%Oligo
    [145 146 7 8 6 4 5 3 2];...%Astro
    [143 144 112 113 114]}; %Micro/Macro
ctypes.Tasic = {[8:16 20:21];[1:7 17:19];24;22;23};

cellgroupname = {'Glut.';'GABA';'Oligo.';'Astro.';'Micro.'};

tkey = ctypes.Tasic;
zkey = ctypes.Zeisel;

clist = lines(length(tkey));

f1 = figure; hold on;
set(f1,'Position',[0 0 1500 290]); hold on; %this to set the size
subplot(1,length(tkey),1,'Parent',f1);
for i = 1:length(tkey)
    curtkey = tkey{i};
    curzkey = zkey{i};
    curBt = tasic.outstruct.Bmeans(:,curtkey);
    Bt = sum(curBt,2);
    maxbt = max(Bt);
    curBz = zeisel.outstruct.Bmeans(:,curzkey);
    Bz = sum(curBz,2);
    maxbz = max(Bz);
    rval(i) = corr(Bt,Bz,'Type','Pearson');
    linr(i) = LinRcalc(Bt,Bz);
    subplot(1,length(tkey),i)
    scatter(Bt,Bz,20,clist(i,:),'filled'); hold on;
    maxlab = max(cat(1,Bt,Bz));
    p = polyfit(Bt,Bz,1);
    x_maxy = (maxlab - p(2))/p(1);
    if x_maxy < maxlab
        plot([0,x_maxy],[p(2),(p(2) + x_maxy*p(1))],'Color',[0 0 0],'LineStyle','--','LineWidth',1.5); hold on;
    else
        plot([0,maxlab],[p(2),(p(2) + maxlab*p(1))],'Color',[0 0 0],'LineStyle','--','LineWidth',1.5); hold on;
    end
    plot([1 max(cat(1,maxbt,maxbz))],[1 max(cat(1,maxbt,maxbz))],'Color',[0 0 0],'LineWidth',1.5); hold on;
    xticks([0 round(max(cat(1,maxbt,maxbz))/4) round(max(cat(1,maxbt,maxbz))/2) 3*(round(max(cat(1,maxbt,maxbz))/4))]);
    yticks([0 round(max(cat(1,maxbt,maxbz))/4) round(max(cat(1,maxbt,maxbz))/2) 3*(round(max(cat(1,maxbt,maxbz))/4))]);
    set(gca,'FontSize',14);
    title({cellgroupname{i};['R = ' num2str(round(rval(i),2)) ', R_C = ' num2str(round(linr(i),2))]},'FontSize',18);
    xlabel('Tasic, et al.','FontSize',18);
    ylabel('Zeisel, et al.','FontSize',18);
    xlim([1 max(cat(1,maxbt,maxbz))]);
    ylim([1 max(cat(1,maxbt,maxbz))]);
end
sgtitle('Mean Cells Per Voxel Across Regions','FontSize',22);

% subplot(3,length(tkey),length(tkey)+1:3*length(tkey));
% bar([rval;linr].','LineWidth',0.0000001); hold on;
% colormap(cbars);
% xticks(1:length(rval));
% xticklabels(cellgroupname);
% yticks([0 0.2 0.4 0.6 0.8 1]);
% set(gca,'FontSize',16);
% ylabel('R & R_C Values','FontSize',20);
% title('Mean Cells Per Voxel: Tasic, et al., Vs. Zeisel, et al.','FontSize',22);
% legend({'R';'R_C'},'Location','northwest');
% resort_ctypes = [1 2 6 3 4 5 7];
% ctypes.Zeisel = [ctypes.Zeisel;{[19:20 43:46 66:71 73:83];...
%     [1 12 15 34 37:42 72 126:127 133:135 138:142 147:149 195:200]}];
% ctypes.Zeisel{1}(8) = [];
% ctypes.Zeisel{3}([6 7]) = [];
% ctypes.Zeisel{5}([1 5]) = [];
% ctypes.Zeisel = ctypes.Zeisel(resort_ctypes);
% cmap = lines(7);
% cmap = cmap(resort_ctypes);

if savenclose
    print('Figure_6f_scatter','-dtiffn');
    close
end
end


