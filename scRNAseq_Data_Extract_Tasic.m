function classstruct = scRNAseq_Data_Extract_Tasic(directory)
% Function that creates a hierarchically organized MATLAB struct object
% that contains all of the cells in the Tasic et al., 2018 and the AIBS,
% 2018 studies grouped by type and region of origin. These data were
% downloaded from the AIBS website:
% https://portal.brain-map.org/atlases-and-data/rnaseq and then extracted
% into a MATLAB file (ReadData.mat).

if nargin < 1
    directory = [cd filesep 'MatFiles'];
end

regions = {'VISp','ALM','LGd'};
load([directory filesep 'ReadData.mat'],'datstruct');
clusterstruct = struct;
classstruct = struct;
for i=1:length(regions)
    structfilename = sprintf('mouse_%s_2018-06-14_samples-columns.mat',regions{i});
    newstructfilename = [directory filesep structfilename];
    load(newstructfilename,'filestruct')
    clusterids = unique(filestruct.cluster);
    highqualind = ones(1,length(clusterids));
    for j=1:length(clusterids)
        if strcmp('Low Q', clusterids{j}(1:5))
            highqualind(j) = 0;
        elseif strcmp('Doub', clusterids{j}(1:4))
            highqualind(j) = 0;
        end
    end
    clusterids = clusterids(logical(highqualind));
    clusterclasslegend = cell(length(clusterids),3);
    clusterclasslegend(:,1) = clusterids;
    for k=1:length(clusterids)
        idcmp = @(x) strcmp(clusterids{k},x);
        idind = cellfun(idcmp, filestruct.cluster);
        clustervarname = strrep(clusterids{k},' ','_');
        clustervarname = strrep(clustervarname,'/','_');
        clusterstruct.(regions{i}).(clustervarname).truename = clusterids{k};
        clusterstruct.(regions{i}).(clustervarname).sampleids = filestruct.sample_name(idind);
        clusterstruct.(regions{i}).(clustervarname).corrs = filestruct.cluster_correlation(idind);
        clusterclasslegend(k,2) = unique(filestruct.subclass(idind));
        clusterclasslegend(k,3) = unique(filestruct.class(idind));
        exonintronind = ismember(datstruct.(regions{i}).Exon.sid, clusterstruct.(regions{i}).(clustervarname).sampleids);
        clusterstruct.(regions{i}).(clustervarname).exonreads = datstruct.(regions{i}).Exon.reads(:,exonintronind);
        clusterstruct.(regions{i}).(clustervarname).intronreads = datstruct.(regions{i}).Intron.reads(:,exonintronind);
        allreads = clusterstruct.(regions{i}).(clustervarname).exonreads + clusterstruct.(regions{i}).(clustervarname).intronreads;
        totalsamplereads = sum(allreads,1);
        normalizedexpr = log2(1000000*allreads./totalsamplereads + 1);
        clusterstruct.(regions{i}).(clustervarname).normalizedexpr = normalizedexpr;
        sids = datstruct.(regions{i}).Exon.sid(exonintronind);
        [test, reorderind] = ismember(sids, clusterstruct.(regions{i}).(clustervarname).sampleids); 
        if all(isnan(clusterstruct.(regions{i}).(clustervarname).corrs))
            clusterstruct.(regions{i}).(clustervarname).centroidexpr = mean(clusterstruct.(regions{i}).(clustervarname).normalizedexpr,2);
        else
            nonanind = ~isnan(clusterstruct.(regions{i}).(clustervarname).corrs);
            corrsnonans = clusterstruct.(regions{i}).(clustervarname).corrs(nonanind);
            exprreorder = clusterstruct.(regions{i}).(clustervarname).normalizedexpr(:,reorderind);
            exprnonans = exprreorder(:,nonanind);
            centroid = exprnonans*corrsnonans/sum(corrsnonans);
            clusterstruct.(regions{i}).(clustervarname).centroidexpr = centroid;
        end
    end
    classids = unique(clusterclasslegend(:,3));
    for l = 1:length(classids)
        classvarname = strrep(classids{l},' ','_');
        classvarname = strrep(classvarname,'/','_');
        classvarname = strrep(classvarname,'-','_');
        classstruct.(regions{i}).(classvarname).truename = classids{l};
        idcmp2 = @(x) strcmp(classids{l},x);
        idind2 = cellfun(idcmp2, clusterclasslegend(:,3));
        subclassids = unique(clusterclasslegend(idind2,2));
        classexprmat = zeros(45768,length(subclassids));
        subclasscellnumber = zeros(1,length(subclassids));
        for m = 1:length(subclassids)
            subclassvarname = strrep(subclassids{m},' ','_');
            subclassvarname = strrep(subclassvarname,'/','_');
            subclassvarname = strrep(subclassvarname,'-','_');
            classstruct.(regions{i}).(classvarname).(subclassvarname).truename = subclassids{m};
            idcmp3 = @(x) strcmp(subclassids{m},x);
            idind3 = cellfun(idcmp3, clusterclasslegend(:,2));
            subclusterids = unique(clusterclasslegend(idind3,1));
            subclassexprmat = zeros(45768,length(subclusterids));
            subclustercellnumber = zeros(1,length(subclusterids));
            for n = 1:length(subclusterids)
                subclustervarname = strrep(subclusterids{n},' ','_');
                subclustervarname = strrep(subclustervarname,'/','_');
                classstruct.(regions{i}).(classvarname).(subclassvarname).(subclustervarname).truename = clusterstruct.(regions{i}).(subclustervarname).truename;
                classstruct.(regions{i}).(classvarname).(subclassvarname).(subclustervarname).sampleids = clusterstruct.(regions{i}).(subclustervarname).sampleids;
                classstruct.(regions{i}).(classvarname).(subclassvarname).(subclustervarname).corrs = clusterstruct.(regions{i}).(subclustervarname).corrs;
                classstruct.(regions{i}).(classvarname).(subclassvarname).(subclustervarname).exonreads = clusterstruct.(regions{i}).(subclustervarname).exonreads;
                classstruct.(regions{i}).(classvarname).(subclassvarname).(subclustervarname).intronreads = clusterstruct.(regions{i}).(subclustervarname).intronreads;
                classstruct.(regions{i}).(classvarname).(subclassvarname).(subclustervarname).normalizedexpr = clusterstruct.(regions{i}).(subclustervarname).normalizedexpr;
                classstruct.(regions{i}).(classvarname).(subclassvarname).(subclustervarname).centroidexpr = clusterstruct.(regions{i}).(subclustervarname).centroidexpr;
                subclustercellnumber(n) = size(classstruct.(regions{i}).(classvarname).(subclassvarname).(subclustervarname).exonreads,2);
                classstruct.(regions{i}).(classvarname).(subclassvarname).(subclustervarname).cellcount = subclustercellnumber(n);
                subclassexprmat(:,n) = classstruct.(regions{i}).(classvarname).(subclassvarname).(subclustervarname).centroidexpr;
            end
            cellnormmeanexpr = subclassexprmat * subclustercellnumber.' / sum(subclustercellnumber);
            classstruct.(regions{i}).(classvarname).(subclassvarname).meanexpr = cellnormmeanexpr;
            subclasscellnumber(m) = sum(subclustercellnumber);
            classexprmat(:,m) = cellnormmeanexpr * subclasscellnumber(m);
            classstruct.(regions{i}).(classvarname).(subclassvarname).cellcount = subclasscellnumber(m);
        end
        classstruct.(regions{i}).(classvarname).meanexpr = sum(classexprmat,2)/sum(subclasscellnumber);
        classstruct.(regions{i}).(classvarname).cellcount = sum(subclasscellnumber);
    end
    classstruct.(regions{i}).EntrezID = datstruct.(regions{i}).Exon.eid;
    clear filestruct
end
end
