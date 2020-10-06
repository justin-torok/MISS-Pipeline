function Figure5_ZeiselCompMaps(outstruct,idx,typeinds,slicelocs,savenclose,directory)

if nargin < 6
    directory = [cd filesep 'MatFiles'];
    if nargin < 5
        savenclose = 0;
        if nargin < 4
            slicelocs = [25,30,47];
        end
    end
end

load([directory filesep 'Zeisel_coronal_geneset.mat'],'Zeisel_gene_names');
load([directory filesep 'Zeisel_Inputs.mat'],'meanexprmat','classkey','regvgene','entrez_names');
load([directory filesep 'PresetInputs.mat'],'GENGDmod','nonzerovox');
load([directory filesep 'input_struct_voxelrender.mat'],'input_struct');
typenames = classkey(typeinds);

unitygenes = ismember(entrez_names,Zeisel_gene_names);
E_corr = regvgene(:,unitygenes);
C_corr = meanexprmat(:,unitygenes);
E_corr = E_corr.'; C_corr = C_corr.';
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
    dex = 1;
    figure('Position',[0 0 1250 375]); subplot(4,5,1); hold on;
    
    for i = 1:2
        
        curmap = zeros(size(GENGDmod));
        switch i
            case 1
                curmap(nonzerovox) = Bcorr(:,typeinds(k));
                ylab = 'Corr.';
                dex = 1;
                threshmax = 0.75;
            case 2
                curmap(nonzerovox) = Bdense(:,typeinds(k));
                ylab = 'MISS';
                dex = 11;
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
            im = interpn(im,2,'spline');
            bw = squeeze(newVoxMap(curloc,:,:));
            im_ = imdilate(bw,se);
            im_ = imerode(im_,se);
            bim_ = interpn(im_,2,'spline');
            bim_(bim_ < 0.5) = 0;
            bim_(bim_ >= 0.5) = 1;
            slice_final = im .* bim_;
            [szy,szx] = size(slice_final);
            ycut = floor(szy*0.9);
            slice_final = slice_final(1:ycut,:);
            bwbounds = bwboundaries(bim_);
            
            subplot(4,5,[dex dex+5]); hold on;
            imagesc(slice_final,[0 threshmax*curmax+eps]); hold on;
            colormap(flipud(pink)); hold on;
            for m = 1:length(bwbounds)
                boundary = bwbounds{m};
                plot(boundary(:,2),boundary(:,1),'k','LineWidth',0.5); hold on;
            end
            set(gca,'xtick',[]);
            set(gca,'ytick',[]);
            set(gca,'Ydir','reverse')
            box on;
            set(gca,'BoxStyle','full');
            
            if j == 1
                ylabel(ylab,'FontSize',14,'FontWeight','bold');
            end
            
            if i == 1
                title(['Slice ' num2str(curloc)],'FontSize',14);
            end
            
            dex = dex + 1;
        
        end
    end
    
    xvar = Bdense(:,typeinds(k));
    yvar = Bcorr(:,typeinds(k)); %yvar(yvar<0) = 0;
%     exclinds = (yvar<0);
%     yvar(exclinds) = []; xvar(exclinds) = [];
    rval = corr(xvar,yvar);
%     rho = corr(xvar,yvar,'Type','Spearman');
    subplot(4,5,[4:5 9:10 14:15 19:20]); hold on;
    scatter(xvar,yvar,40,'ko','filled'); hold on;
    bfl = lsline;
    bfl.Color = [0.85 0.33 0.1];
    bfl.LineWidth = 2;
    xlabel('MISS Counts','FontSize',14,'FontWeight','bold');
    yyaxis right
    ylabel('Correlation','FontSize',14,'FontWeight','bold');
    yticks([]);
    yyaxis left
    yticks([0 round(max(yvar)/2,2) round(max(yvar),2)]);
    xticks([0 floor(max(xvar)/3) floor(2*(max(xvar)/3)) floor(max(xvar))]);
%     title(['R = ' num2str(round(rval,2)) ', \rho = ' num2str(round(rho))],'FontSize',14);
    title(['R = ' num2str(round(rval,2))],'FontSize',14);
    set(gca,'FontSize',14);
    
    sgtitle(typenames{k},'FontSize',18);
    
    if savenclose
        print(['Zeisel_MISS_Comp_' typenames{k}],'-dtiffn');
        close
    end
    clear slice_raw slice_final
    
end


