function Figure_6bcde_zeisel_typemaps(outstruct,idx,typeinds,slicelocs,savenclose,directory)

if nargin < 6
    directory = [cd filesep 'MatFiles'];
    if nargin < 5
        savenclose = 0;
        if nargin < 4
            slicelocs = {[35 45 60];[35 40 65];[40 45 65];[35 45 60];[30 50 65]};
            if nargin < 3
                typeinds = [153 27 83 92 117];
            end
        end
    end
end

load([directory filesep 'Zeisel_coronal_geneset.mat'],'Zeisel_gene_names');
load([directory filesep 'Zeisel_Inputs.mat'],'genevct','classkey',...
    'voxvgene','GENGDmod','nonzerovox','gene_names');
load([directory filesep 'input_struct_voxelrender.mat'],'input_struct');
typenames = classkey(typeinds);

unitygenes = ismember(gene_names,Zeisel_gene_names);
E_corr = voxvgene(:,unitygenes).';
C_corr = genevct(unitygenes,:);
% E_corr = E_corr.'; C_corr = C_corr.';
nvox = size(E_corr,2);
logE = log2(E_corr + 1); % Zeisel et al initial log2 normalization
mn_logE = mean(logE,2); std_logE = std(logE,[],2);
Z_logE = (logE - repmat(mn_logE,1,nvox)) ./ repmat(std_logE,1,nvox);
Enans = isnan(Z_logE);
ntypes = size(C_corr,2);
mn_C = mean(C_corr,2); 
Z_C = (C_corr - repmat(mn_C,1,ntypes)) ./ repmat(mn_C,1,ntypes);
Cnans = isnan(Z_C);
allnans = logical(Cnans(:,1) + Enans(:,1)); 
Bcorr = corr(Z_C(~allnans,:),Z_logE(~allnans,:));
Bcorr = Bcorr.';

Bdense = outstruct(idx).corrB;

% maplocs = slicelocs*2 - 1;
newVoxMap = input_struct.brain_atlas;
newVoxMap = double(logical(newVoxMap));

for k = 1:length(typeinds)
    
    maplocs = slicelocs{k};
    maplocs = maplocs*2 - 1;
%     dex = 1;
    figure('Position',[0 0 400 1200]); subplot(10,2,1); hold on;
    
    for i = 1:2
        
        curmap = zeros(size(GENGDmod));
        switch i
            case 1
                curmap(nonzerovox) = Bcorr(:,typeinds(k));
                tlab = 'Corr.';
                dex = 1;
                threshmax = 0.75;
            case 2
                curmap(nonzerovox) = Bdense(:,typeinds(k));
                tlab = 'MISS';
                dex = 2;
                threshmax = 0.5;
        end
        curmap = imresize3(curmap,[133 81 115]);
        curmap(curmap<0) = 0;
        curmax = max(max(max(curmap)));
    
        for j = 1:length(maplocs)
            
            curloc = maplocs(j);
            slice_raw = squeeze(curmap(curloc,:,:));
            im = slice_raw;
            se = offsetstrel('ball',3,1,4);
            im = imdilate(im,se);
            im = imerode(im,se);
            im = interpn(im,1,'spline');
            bw = squeeze(newVoxMap(curloc,:,:));
            im_ = imdilate(bw,se);
            im_ = imerode(im_,se);
            bim_ = interpn(im_,1,'spline');
            bim_(bim_ < 0.5) = 0;
            bim_(bim_ >= 0.5) = 1;
            slice_final = im .* bim_;
%             [szy,szx] = size(slice_final);
%             ycut = floor(szy*0.9);
%             slice_final = slice_final(1:ycut,:);
            bwbounds = bwboundaries(bim_);
            plotmaxes = zeros(length(bwbounds),2); plotmins = plotmaxes;
            for m = 1:length(bwbounds)
                plotmaxes(m,:) = max(bwbounds{m});
                plotmins(m,:) = min(bwbounds{m});
            end
            plotmaxes = max(plotmaxes);
            plotmins = min(plotmins);
            subplot(10,2,[dex dex+2]); hold on;
            imagesc(slice_final,[0 threshmax*curmax+eps]); hold on;
            colormap(flipud(pink)); hold on;
            for m = 1:length(bwbounds)
                boundary = bwbounds{m};
                plot(boundary(:,2),boundary(:,1),'k','LineWidth',0.5); hold on;
            end
            ylim([plotmins(1)-5, plotmaxes(1)+5]); xlim([plotmins(2)-5, plotmaxes(2)+5]);
            set(gca,'XAxisLocation','origin')
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            set(gca,'Ydir','reverse')
            box on;
            set(gca,'BoxStyle','full');
            if i == 1
                ylabel(['Slice ' num2str((curloc+1)/2)],'FontSize',18,'FontWeight','bold');
            end
            
            if j == 1
                title(tlab,'FontSize',18);
            end
            
            dex = dex + 4;
        
        end
    end
    
    xvar = Bdense(:,typeinds(k));
    yvar = Bcorr(:,typeinds(k)); %yvar(yvar<0) = 0;
%     exclinds = (yvar<0);
%     yvar(exclinds) = []; xvar(exclinds) = [];
    rval = corr(xvar,yvar);
%     rho = corr(xvar,yvar,'Type','Spearman');
    subplot(10,2,13:20); hold on;
    scatter(xvar,yvar,40,'ko','filled'); hold on;
    bfl = lsline;
    bfl.Color = [0.85 0.33 0.1];
    bfl.LineWidth = 2;
    xlabel('MISS Counts','FontSize',18,'FontWeight','bold');
    yyaxis right
    ylabel('Correlation','FontSize',18,'FontWeight','bold');
    yticks([]);
    yyaxis left
    yticks([0 round(max(yvar)/2,2) round(max(yvar),2)]);
    xticks([0 floor(max(xvar)/3) floor(2*(max(xvar)/3)) floor(max(xvar))]);
%     title(['R = ' num2str(round(rval,2)) ', \rho = ' num2str(round(rho))],'FontSize',14);
%     title(['R = ' num2str(round(rval,2))],'FontSize',14);
    txt = ['R = ' num2str(round(rval,2))];
    text(0.1*max(xvar),1.15*max(yvar),txt,'FontSize',18,'FontWeight','bold');
    set(gca,'FontSize',14);
    
    sgtitle(typenames{k},'FontSize',18);
    
    if savenclose
        print(['Zeisel_MISS_Comp_' typenames{k}],'-dtiffn');
        close
    end
    clear slice_raw slice_final
    
end


