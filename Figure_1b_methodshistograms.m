function Figure_1b_methodshistograms(method,ngen_param,typeinds,savenclose,directory)
%all non-ontology vars necessary for generating input matrices

% strrepper = @(x) strrep(x, '_', ' ');
if nargin < 5
    directory = [cd filesep 'MatFiles'];
    if nargin < 4
        savenclose = 0;
    end
end
load([directory filesep 'PresetInputs.mat'],'classkey','meanexprmat','entrez_names','regvgene');

%pre-diff exp analysis necessary normalizations
ctmean = mean(meanexprmat,1);
ctmean = repmat(ctmean,size(meanexprmat,1),1);
ctnorm = meanexprmat ./ ctmean;
ctnorm(ctnorm<(0.1*std(nonzeros(ctnorm)))) = 0;

genesum = sum(ctnorm,2);
genesum = repmat(genesum,1,size(meanexprmat,2));
genenorm = ctnorm ./ genesum;
isnans = isnan(genenorm(:,1));
genenorm(isnans,:) = 0;

%finding diff exp genes, getting names, weeding out duplicates according to
%different methods. I'm keeping the code longer than strictly necessary to
%improve readability
cutoffs = zeros(1,25);
if strcmp(method,'COLAMD')
    thresh = 0.06; %making matrix x% sparse before colamd
    normthresh = genenorm;
    normthresh(normthresh<thresh) = 0;
    normthresh = normthresh.';
    CC = colamd(normthresh);
    cainds = CC(1:ngen_param);
    entropyscore = zeros(1,size(genenorm,1));
    for i = 1:size(genenorm,1)
        curtype = genenorm(i,:);
        curtype(curtype == 0) = 1; % necessary for entropy calculation
        entropyscore(i) = -dot(log(curtype),curtype);
    end
    genset_names = entrez_names(cainds);
elseif strcmp(method,'DBSCAN')
    genset_names = 'a';
    epsilon = ngen_param; % distance parameter for DBSCAN
    for i = 1:size(genenorm,2)
        curtype = genenorm(:,i);
        dbinds = dbscan(curtype,epsilon,30);
        genenames{i} = entrez_names(dbinds == -1);
        genenames{i} = sort(genenames{i});
        genset_names = [genset_names;genenames{i}];
        cutoffs(i) = min(curtype(dbinds == -1));
        DBinds(:,i) = dbinds;
    end
    genset_names(1) = [];
    genset_names = unique(genset_names);
elseif strcmp(method,'Entropy')
    entthresh = ngen_param; % entropy score threshold
    entropyscore = zeros(1,size(genenorm,1));
    for i = 1:size(genenorm,1)
        curtype = genenorm(i,:);
        curtype(curtype == 0) = 1; % necessary for entropy calculation
        entropyscore(i) = -dot(log(curtype),curtype);
    end
    entinds = find(entropyscore < entthresh);
    zeroinds = find(isnans);
    cutoff = max(entropyscore(entinds));
    genset_names = entrez_names(setdiff(entinds,zeroinds));
elseif strcmp(method,'mRMR')
    mrmr_n = ngen_param;
    entropyscore = zeros(1,size(genenorm,1));
    for i = 1:size(genenorm,1)
        curtype = genenorm(i,:);
        curtype(curtype == 0) = 1; % necessary for entropy calculation
        entropyscore(i) = -dot(log(curtype),curtype);
    end
    mrmrinds = mRMR_Selector(ctnorm,mrmr_n,'Quo');
    genset_names = entrez_names(mrmrinds);
elseif strcmp(method,'MRx3')
    mrmr_n = ngen_param;
    lambda = 150;
    entropyscore = zeros(1,size(genenorm,1));
    for i = 1:size(genenorm,1)
        curtype = genenorm(i,:);
        curtype(curtype == 0) = 1; % necessary for entropy calculation
        entropyscore(i) = -dot(log(curtype),curtype);
    end
    [mrmrinds] = MRx3_Selector(meanexprmat,regvgene,mrmr_n,lambda);
    mrmrinds = sort(mrmrinds);
    genset_names = entrez_names(mrmrinds);
end

%generating indices for diff exp genes, choosing in spatexp and meanexprmat
gendex = ismember(entrez_names,genset_names);
notgendex = ~gendex; 
notgendex = logical(notgendex - isnans);

if strcmp(method,'Entropy')
    figure; hold on;
    h1 = histogram(entropyscore(notgendex));
    h1.BinWidth = 0.075;
    h1.FaceColor = [1 1 1];
    h2 = histogram(entropyscore(gendex));
    h2.BinWidth = 0.075;
    h2.FaceColor = [0 0 0];
    plot(cutoff, 300, 'r*', 'MarkerSize', 10);
    set(gca,'YLim',[0,300],'XLim',[0,max(entropyscore)]);
    set(gca,'XTick',[],'YTick',[],'FontSize',40);
    if savenclose
        saveas(gcf,[method '_genediscrim_hist'],'tiffn');
        close
    end
elseif strcmp(method,'COLAMD')
    figure; hold on;
    h1 = histogram(entropyscore(notgendex));
    h1.BinWidth = 0.05;
    h1.FaceColor = [1 1 1];
    h2 = histogram(entropyscore(gendex));
    h2.BinWidth = 0.05;
    h2.FaceColor = [0 0 0];
    set(gca,'YLim',[0,300],'XLim',[0,max(entropyscore)]);
    set(gca,'XTick',[],'YTick',[],'FontSize',40);
    if savenclose
        saveas(gcf,[method '_genediscrim_hist'],'tiffn');
        close
    end
elseif strcmp(method,'mRMR') || strcmp(method,'MRx3')
    figure; hold on;
    h1 = histogram(entropyscore(notgendex));
    h1.BinWidth = 0.075;
    h1.FaceColor = [1 1 1];
    h2 = histogram(entropyscore(gendex));
    h2.BinWidth = 0.075;
    h2.FaceColor = [0 0 0];
    set(gca,'YLim',[0,300],'XLim',[0,max(entropyscore)]);
    set(gca,'XTick',[],'YTick',[],'FontSize',40);
    if savenclose
        saveas(gcf,[method '_genediscrim_hist'],'tiffn');
        close
    end
elseif strcmp(method,'DBSCAN')
    for i = typeinds
%         type = strrepper(classkey{i});
        figure; hold on;
        curtype = genenorm(:,i);
        h2 = histogram(curtype(DBinds(:,i) ~= -1));
        h2.BinWidth = 0.0225;  
        h2.FaceColor = [1 1 1];
        h1 = histogram(curtype(DBinds(:,i) == -1));
        h1.BinWidth = 0.0225;
        h1.FaceColor = [0 0 0];
        plot(cutoffs(i), 150, 'r*', 'MarkerSize', 10);
        set(gca,'YLim',[0,150],'XLim',[0,1]);
        set(gca,'XTick',[],'YTick',[],'FontSize',40);
        if savenclose
            saveas(gcf,[method '_genediscrim_hist_' classkey{2,i}],'tiffn');
            close
        end
    end
end


% 
end