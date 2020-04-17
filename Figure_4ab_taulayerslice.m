function Figure_4ab_taulayerslice(outstruct,idx,slicelocs,savenclose,directory)

if nargin < 5
    directory = [cd filesep 'MatFiles'];
    if nargin < 4
        savenclose = 0;
        if nargin < 3
            slicelocs = [25,31,36];
        end
    end
end
load([directory filesep 'tau_calc_dependencies.mat'],'GENGDmod',...
    'nonzerovox','structList','structIndex','neomask3d','neobounds3d');

ranks = [1 2 3 3 4 4 4];
cell_names = {'L2n3','L4','L5IT','L5PT','L6CT','L6IT','L6b'};
cell_inds = 9:15;
titlecellnames = {'L2/3-IT', 'L4','L5-IT', 'L5-PT', 'L6-CT','L6-IT','L6b'};
taustruct = TauCalc(outstruct,idx,cell_names,cell_inds,ranks,directory);
B = outstruct(idx).Bvals;
ngen = outstruct(idx).nGen;
if ismember('sumfit',fieldnames(outstruct(idx)))
    sumfit = outstruct(idx).sumfit;
else
    sumfit = outstruct(idx).LinR.micro(3) + outstruct(idx).LinR.pv(3)+...
        outstruct(idx).LinR.sst(3)+outstruct(idx).LinR.vip(3) + taustruct.tau;
end
if ismember('method',fieldnames(outstruct(idx)))
    method = [outstruct(idx).method ' '];
else
    method = '';
end

newVoxMap = zeros(size(GENGDmod));
pervoxglut = struct;
for i = 1:length(cell_names)
    pervoxglut.(cell_names{i}) = newVoxMap;
end
newVoxMap(nonzerovox) = 1;
glutnames = cell_names;

for i = 1:length(structList)
    [~,voxinds] = ismember(structIndex{i},nonzerovox);
    allvox = nonzerovox(voxinds);
    curvox = B(voxinds,cell_inds);
    for k = 1:length(glutnames)
        curglut = squeeze(curvox(:,k));
        pervoxglut.(glutnames{k})(allvox) = curglut;
        clear curglut
    end
    clear allvox curvox voxinds
end

[d1,~,~] = find(neomask3d);
d1inds = unique(d1);
maplocs = d1inds;
for k = 1:length(glutnames)
    curmap = pervoxglut.(glutnames{k});
    for j = 10:length(maplocs)
        curloc = maplocs(j);
        ctxsurf = squeeze(neobounds3d(curloc,:,:));
        curneomask = squeeze(neomask3d(curloc,:,:));
        [ctxsurfx,ctxsurfy] = find(ctxsurf);
        curslc = squeeze(curmap(curloc,:,:));
        mnzs = min(nonzeros(curslc(:)));
        if isempty(mnzs)
            curslc(curslc==0) = 0;
        else
            curslc(curslc==0) = mnzs;
        end
        grayslc = mat2gray(curslc);
        binslc_whole = logical(grayslc);
        paleoinds = find(~curneomask);
        binslc = binslc_whole;
        binslc(paleoinds) = 0;
        skeleton = bwskel(binslc);
        [skelx,skely] = find(skeleton);
        for s = 1:length(find(skeleton))
            curx = skelx(s);
            cury = skely(s);
            dists(s) = min(sqrt((ctxsurfx-curx).^2 + (ctxsurfy-cury).^2));
        end
        if isempty(skelx) || isempty(curx) || isempty(dists)
            mndists(j-9,k) = -1;
        else
            mndists(j-9,k) = mean(dists(:));
        end
        clear dists cury curx
    end
end

lranks = 0;
trueranks = 0;
for j = 1:size(mndists,1)
    curmns = mndists(j,:);
    nosig = (curmns==-1);
    n_nosig(j) = length(find(nosig));
    curmns(nosig) = [];
    [~,rankinds] = sort(curmns,'ascend');
    curranks = ranks;
    curranks(nosig) = [];
    sliceranks = curranks;
    trueranks = [trueranks curranks];
    curranks = curranks(rankinds);
    lranks = [lranks curranks];
    if isempty(curranks)
        tau_per = 0;
    else
        [tau_per,~] = corr(curranks.', sliceranks.', 'type', 'Kendall');
    end
    tau_slice(j) = (1-(n_nosig(j)/7))*tau_per;
    sqerrvec(j) = (1 - tau_per);
end

f1 = figure; hold on;
set(f1,'Position',[0 0 350 1200]); hold on; %this to set the size
s1 = subplot(19,3,1,'Parent',f1);
subplot(19,3,1:9); hold on;
plot(maplocs(10):maplocs(end),tau_slice,'LineWidth',2.5); hold on;
set(gca,'FontSize',12);
xlabel('Slice No.','FontSize',14);
ylabel('\tau-Value','FontSize',14);
title([method '\tau-Value Across Slices'],'FontSize',14);
txt = {['\tau-Val = ' sprintf('%.2f',taustruct.tau)];...
    sprintf('Sum Fit = %.2f',sumfit);['{\it n}_G = ' sprintf('%d',ngen)]};
ylim([-0.5 1]);
text(21,0,txt);

redcellnames = glutnames;
slice_name = {'Rostral', 'Middle', 'Caudal'};
maplocs = slicelocs;

index = 16; index2 = 19;
for k = 1:length(redcellnames)
    curmap = pervoxglut.(redcellnames{k});
    for j = 1:length(maplocs)
        curloc = maplocs(j);
        slice_raw = squeeze(curmap(curloc,:,:));
        curneo = squeeze(neomask3d(curloc,:,:));
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
        curmin = min(min(nonzeros(slice_final)));
        bw = squeeze(newVoxMap(curloc,:,:));
        im_ = imdilate(bw, se);
        im_ = imerode(im_, se);
        bim_ = interpn(im_,2,'linear');
        bim_(bim_ < 0.67) = 0;
        bim_(bim_ >= 0.67) = 1;
        biminds = logical(bim_(:));
        slice_final(biminds) = slice_final(biminds) + curmax/30;
        bwbounds = bwboundaries(bim_);
        figure(f1); subplot(19,3,[index index2]); hold on;
        imagesc(slice_final,[0 curmax/2]); hold on;
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
        if index == 16
            ylabel(titlecellnames{k},'FontSize',14);
            title(slice_name{j},'FontSize',14);
        elseif index == 17 || index == 18
            title(slice_name{j},'FontSize',14);
        elseif index == 22 || index == 28 || index == 34
            ylabel(titlecellnames{k},'FontSize',14);
        elseif index == 40 || index == 46 || index == 52
            ylabel(titlecellnames{k},'FontSize',14);
        end
        index = index + 1;
        index2 = index2 + 1;
        clear slice_raw slice_final
    end
    index = index + 3;
    index2 = index2 + 3;
end

if savenclose
    print('Figure_4ab_layertype','-dtiff');
    close
end
end