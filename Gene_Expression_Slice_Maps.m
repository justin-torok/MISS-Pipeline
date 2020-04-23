function Gene_Expression_Slice_Maps(gene_names,slicelocs,savenclose,directory)

% 1. Load raw data and dependencies for mapping; define defaults
if nargin < 4
    directory = [cd filesep 'MatFiles'];
    if nargin < 3
        savenclose = 0;
        if nargin < 2
            slicelocs = 34;
        end
    end
end
load([directory filesep 'PresetInputs.mat'],'regvgene','entrez_names','GENGDmod','nonzerovox');
% 2. Define density matrix and ensure gene names are valid
test = ismember(gene_names,entrez_names);
if ~all(test)
    invalidgenes = gene_names(~test);
    invalidgenesstr = invalidgenes{1};
    i = 1;
    while i < length(invalidgenes)
        i = i + 1;
        invalidgenesstr = [invalidgenesstr ', ' invalidgenes{i}];
    end
    error('Invalid gene names: %s', invalidgenesstr)
end
geneinds = find(ismember(entrez_names,gene_names));
D = regvgene;

% 3. Map the results

maplocs = slicelocs;
newVoxMap = zeros(size(GENGDmod));
newVoxMap(nonzerovox) = 1;

for k = 1:length(gene_names)
    curmap = zeros(size(GENGDmod));
    curmap(nonzerovox) = D(:,geneinds(k));
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
        title(sprintf('%s Density Map',gene_names{k}));
        set(gca,'FontSize',24);
        if savenclose
            print(sprintf('GeneExprDensityMap_%s_Slice_%d',gene_names{k},maplocs(j)),'-dtiffn');
            close
        end
        clear slice_raw slice_final
    end
end

end

