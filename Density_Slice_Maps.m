function Density_Slice_Maps(outstruct,idx,typeinds,slicelocs,savenclose,study,directory)

% 1. Load raw data and dependencies for mapping; define defaults
if nargin < 7
    directory = [cd filesep 'MatFiles'];
    if nargin < 6
        study = 'tasic';
        if nargin < 5
            savenclose = 0;
            if nargin < 4
                slicelocs = [25,30,47];
            end
        end
    end
end
if strcmp(study,'tasic')
    load([directory filesep 'PresetInputs.mat'],'classkey','GENGDmod','nonzerovox');
elseif strcmp(study,'zeisel')
    load([directory filesep 'PresetInputs.mat'],'GENGDmod','nonzerovox');
    load([directory filesep 'Zeisel_extract.mat'],'classkey');
end
typenames = classkey(typeinds);

% 2. Define densities
D = outstruct(idx).corrB;

% 3. Map the results

maplocs = slicelocs;
newVoxMap = zeros(size(GENGDmod));
newVoxMap(nonzerovox) = 1;

for k = 1:length(typeinds)
    curmap = zeros(size(GENGDmod));
    curmap(nonzerovox) = D(:,typeinds(k));
    curmax = max(max(max(curmap)));
    for j = 1:length(maplocs)
        curloc = maplocs(j);
        slice_raw = squeeze(curmap(curloc,:,:));
        im = slice_raw;
        se = strel('disk',1,4);
        im = imdilate(im, se);
        im = imerode(im, se);
        bim = double(logical(im));
        im = interpn(im,3,'spline');
        bim = interpn(bim,3,'spline');
        bim(bim < 0.67) = 0;
        bim(bim >= 0.67) = 1;
        slice_final = im .* bim;
%         curmax = max(max(slice_final));
        bw = squeeze(newVoxMap(curloc,:,:));
        im_ = imdilate(bw, se);
        im_ = imerode(im_, se);
        bim_ = interpn(im_,3,'spline');
        bim_(bim_ < 0.67) = 0;
        bim_(bim_ >= 0.67) = 1;
        biminds = logical(bim_(:));
%         slice_final(biminds) = slice_final(biminds) + curmax/30;
        bwbounds = bwboundaries(bim_);
        
        f1 = figure;
        set(f1,'Position',[0 0 600 500]); %this to set the size
        imagesc(slice_final,[0 0.6*curmax+eps]); hold on;
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
        title(sprintf('%s Density Map, nG = %d',typenames{k},outstruct(idx).nGen));
        set(gca,'FontSize',24);
        if savenclose
            print(sprintf('DensityMap_%s_Slice_%d',typenames{k},maplocs(j)),'-dtiffn');
            close
        end
        clear slice_raw slice_final
    end
end

end

