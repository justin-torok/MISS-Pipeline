function Figure_2d_correlationmaps(typeinds,method,ngen_param,lambda,slicelocs,savenclose,study,directory)

% 1. Load raw data and dependencies for mapping; define defaults; create
% E_red, C_red
if nargin < 8
    directory = [cd filesep 'MatFiles'];
    if nargin < 7
        study = 'tasic';
        if nargin < 6
            savenclose = 0;
            if nargin < 5
                slicelocs = [25,30,47];
                if nargin < 4
                    lambda = 150;
                end
            end
        end
    end
end
if strcmp(study,'tasic')
    load([directory filesep 'PresetInputs.mat'],'meanexprmat','regvgene',...
        'entrez_names','classkey','GENGDmod','nonzerovox');
    [E,C] = GeneSelector(meanexprmat,regvgene,entrez_names,ngen_param,lambda,method);
    typenames = classkey(typeinds);
elseif strcmp(study,'zeisel')
    load([directory filesep 'Zeisel_extract.mat'],'meanexprmat','entrez_names','classkey');
    regvgene = ISH_Data_Extract_Zeisel(directory);
    meanexprmat = meanexprmat.';
    [E,C] = GeneSelector(meanexprmat,regvgene,entrez_names,ngen_param,lambda,method);
    typenames = classkey(typeinds);
    load([directory filesep 'PresetInputs.mat'],'GENGDmod','nonzerovox');
end

% 2. Renormalize the ISH data according to Zeisel et al, 2018
% This does not include a pre-filtering step; we assume that our MRx3
% procedure removes genes that are too low-signal for our purposes from
% consideration, as that incorporates a per-gene penalty for poor cell
% density inference.
nvox = size(E,2);
logE = log2(E + 1); % Zeisel et al initial log2 normalization
mn_logE = mean(logE,2); std_logE = std(logE,[],2);
Z_logE = (logE - repmat(mn_logE,1,nvox)) ./ repmat(std_logE,1,nvox);
Enans = isnan(Z_logE);

% 3. Renormalize the scRNAseq data
% Note that C has already been log2-transformed
ntypes = size(C,2);
mn_C = mean(C,2); std_C = std(C,[],2);
Z_C = (C - repmat(mn_C,1,ntypes)) ./ repmat(mn_C,1,ntypes);
Cnans = isnan(Z_C);

% 4. Determine the voxels x types cross-correlation matrix, Bcorr
allnans = logical(Cnans(:,1) + Enans(:,1)); 
Bcorr = corr(Z_C(~allnans,:),Z_logE(~allnans,:));
Bcorr = Bcorr.';

% 5. Map the results

maplocs = slicelocs;
newVoxMap = zeros(size(GENGDmod));
newVoxMap(nonzerovox) = 1;

for k = 1:length(typeinds)
    curmap = zeros(size(GENGDmod));
    curmap(nonzerovox) = Bcorr(:,typeinds(k));
    for j = 1:length(maplocs)
        curloc = maplocs(j);
        slice_raw = squeeze(curmap(curloc,:,:));
        im = slice_raw;
        se = strel('diamond',1);
        im = imdilate(im, se);
        im = imerode(im, se);
        bim = double(logical(im));
        im = interpn(im,2,'linear');
        bim = interpn(bim,2,'linear');
        bim(bim < 0.67) = 0;
        bim(bim >= 0.67) = 1;
        slice_final = im .* bim;
        curmax = max(max(slice_final));
        bw = squeeze(newVoxMap(curloc,:,:));
        im_ = imdilate(bw, se);
        im_ = imerode(im_, se);
        bim_ = interpn(im_,2,'linear');
        bim_(bim_ < 0.67) = 0;
        bim_(bim_ >= 0.67) = 1;
        biminds = logical(bim_(:));
        slice_final(biminds) = slice_final(biminds) + curmax/30;
        bwbounds = bwboundaries(bim_);
        
        f1 = figure;
        set(f1,'Position',[0 0 600 500]); %this to set the size
        imagesc(slice_final,[0 0.85*curmax]); hold on;
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
        title(sprintf('%s Z-score Correlation, nG = %d',typenames{k},size(C,1)));
        set(gca,'FontSize',24);
        if savenclose
            print(sprintf('CorrelationMap_%s_Slice_%d',typenames{k},maplocs(j)),'-dtiffn');
            close
        end
        clear slice_raw slice_final
    end
end

end

