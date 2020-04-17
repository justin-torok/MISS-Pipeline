function brainframe(input_struct)

%Voxel or region resolution flag
voxUreg = input_struct.voxUreg;

%Base brain atlas prep
brainat = input_struct.brain_atlas;
CCFinds = brainat;
CCFbin = brainat;
CCFbin(CCFbin > 0) = 1;
isomap = CCFbin;
isomap = smooth3(isomap,'box',5);
isobin = isomap;
isobin(isobin > 0) = 1;
CCFinds = CCFinds .* isobin;
surf1 = isosurface(isomap);
p1 = patch(surf1);
% v = get(p1);
% w = isonormals(isomap,p1);
if strcmp(input_struct.bgcolor,'w')
    set(p1,'FaceColor',[0.9 0.85 0.85],'EdgeColor','none','FaceAlpha',0.05);
elseif strcmp(input_struct.bgcolor,'k')% set the color, mesh and transparency level of the surface
    set(p1,'FaceColor',[1 1 1],'EdgeColor','none','FaceAlpha',0.05);
else
    set(p1,'FaceColor',[0.95 0.9 0.9],'EdgeColor','none','FaceAlpha',0.0625);
end
daspect([1,1,1])
view(3); axis tight
camlight; lighting gouraud
hold on

if voxUreg
    %Per-region rendering
    %Setting up data, per region
    groupid = input_struct.region_groups;
    idlist = unique(groupid);
    is0 = (idlist==0);
    idlist(is0) = [];
    scalefac = input_struct.xfac;
    cmap = input_struct.cmap;
    scalevec = input_struct.data;
    normvec = scalevec / mean(scalevec);
    normvec = normvec / scalefac;
    ptsz = input_struct.pointsize;
    sphsz = input_struct.size_constant;
    
    %Grouping regions for plotting
    for i = 1:length(idlist)
        curid = idlist(i);
        reginds = find(groupid==curid);
        centers = [0 0 0];
        normvals = normvec(reginds);
        
        for k = 1:length(reginds)
            if scalevec(k) > 0
                curreg = reginds(k);
                pcinds = find(CCFinds==curreg);
                [x,y,z] = ind2sub(size(CCFbin),pcinds);
                centroid = [mean(y) mean(x) mean(z)];
                %Finding region centers, placing a sphere around them
                if input_struct.sphere
                    [sphx,sphz,sphy] = sphere(ceil(scalefac*scalevec(k)));
                    centroid = repmat(centroid,ceil(scalefac*scalevec(k)+1)^2,1) +...
                        (sphsz*scalevec(k))*[sphx(:) sphy(:) sphz(:)];
                %Setting up diffuse random point clouds per region
                elseif ~input_struct.sphere && ~input_struct.centered(1)
                    rng(k);
                    randinds = randi(length(pcinds),ceil(normvals(k)*length(pcinds)^(1/3)),1);
%                     randinds = randi(length(pcinds),ceil(normvals(k)),1);
                    chosevox = [y(randinds) x(randinds) z(randinds)];
                    centroid = chosevox + rand(size(chosevox,1),3);
                %Creating random point clouds biased towards region centers
                elseif ~input_struct.sphere && input_struct.centered(1)
                    if length(input_struct.centered) > 1
                        centstrength = input_struct.centered(2);
                    else
                        centstrength = 0.5;
                    end
                    voxdists = ((y-centroid(1)).^2 + (x-centroid(2)).^2 + (z-centroid(3)).^2).^(centstrength);
                    dieinds = WeightedDie(voxdists);
%                     npts = ceil(normvals(k));
                    if mean(nonzeros(scalevec)) == 1
%                         npts = ceil(normvals(k)*length(pcinds));
                          npts = ceil(normvals(k));
                    else
                        npts = ceil(normvals(k)*length(pcinds)^(1/3));
                    end
                    rng(k);
                    choseinds = randi(length(dieinds),npts,1);
                    choseinds = dieinds(choseinds);
                    choseinds = pcinds(choseinds);
                    [x_chose,y_chose,z_chose] = ind2sub(size(CCFbin),choseinds);
                    chosevox = [y_chose x_chose z_chose];
                    centroid = chosevox + rand(size(chosevox,1),3); %.* multipliers;
                end
                centers = [centers;centroid];
            end
        end
        
        %Plotting regions at center, rendered as sphere
        centers(1,:) = [];
        dists = sqrt(centers(:,1).^2 + centers(:,2).^2 + centers(:,3).^2);
        intensities = linspace(0.75,1,size(centers,1))';
        [~,sortinds] = sort(dists);          
        intensities = intensities(sortinds);
        ptcloud = pointCloud(centers,'Color',...
            repmat(cmap(i,:),size(centers,1),1),...
            'Intensity',intensities);
        pcshow(ptcloud,'MarkerSize',ptsz); hold on;
    end
else
    %Per-voxel rendering
    %Data mapping and prep, per voxel
    interpmap = input_struct.data;
    testvals = interpmap(:);
    sumvals = sum(testvals);
    sortvals = sort(testvals,'descend');
    cumsumvals = cumsum(sortvals);
    fracvals = cumsumvals / sumvals;
    interpmap = CCFbin .* interpmap;
    cmap = input_struct.cmap;
    interpmap_f = smooth3(interpmap,'box',[7 1 1]);
    interpmap_f = CCFbin .* interpmap_f;
    
    %Setting up binning and misc for visualization
    nbin = input_struct.nbin;
    rangei = linspace(0,0.7,nbin+1);
    rangei = fliplr(rangei);
    xfac = input_struct.xfac;
    ptsz = input_struct.pointsize;
    
    %Looping through bins and visualizing point clouds
    for i = 1:nbin
        bounds = sortvals(fracvals<rangei(i) & fracvals>=rangei(i+1));
        bound1 = bounds(1); bound2 = bounds(end);
        [x,y,z] = ind2sub(size(interpmap_f),find(interpmap_f>=bound2 & interpmap_f<bound1));
        xyzmap = [y x z];
        ifac = xfac * i;
        xyz_jitter = repmat(xyzmap,floor(ifac)+1,1) + [zeros(size(xyzmap));rand(size(xyzmap,1)*floor(ifac),size(xyzmap,2))];
        ptcloud = pointCloud(xyz_jitter,'Color',repmat(cmap(i,:),size(xyz_jitter,1),1),...
            'Intensity',repmat(0.1*i,size(xyz_jitter,1),1));
        pcshow(ptcloud,'MarkerSize',ptsz); hold on;
        clear x y z xyzmap ptcloud xyz_jitter
    end
end

%Changing axis properties for visualization
set(gcf,'color',input_struct.bgcolor);
set(gca,'color',input_struct.bgcolor);
ax = gca;
set(ax,'XColor','none','YColor','none','ZColor','none');

%Saving and closing or opening .fig GUI
savenclose = input_struct.savenclose;
imglab = input_struct.img_labels;
imgtype = input_struct.img_format;
if savenclose
    view([0, 0, 1]);
    ax = gca;
    set(ax,'XColor','none','YColor','none','ZColor','none');
    set(ax,'XTick',[],'YTick',[],'ZTick',[]);
    saveas(gcf,[imglab '_sagittal'],imgtype);
    view([-1, 0, 0]);
    ax = gca;
    set(ax,'XColor','none','YColor','none','ZColor','none');
    set(ax,'XTick',[],'YTick',[],'ZTick',[]);
    saveas(gcf,[imglab '_axial'],imgtype);
    view([0, -1, 0]);
    ax = gca;
    set(ax,'XColor','none','YColor','none','ZColor','none');
    set(ax,'XTick',[],'YTick',[],'ZTick',[]);
    saveas(gcf,[imglab '_coronal'],imgtype);
    close
    clear isomap
end

end

function [weighted_die] = WeightedDie(distarr)
% This function creates a weighted die of indices based on an input of
% distances, where the "weight" is inversely proportional to the distance.

% 1. Create a "similarity" capped at 100 with a minimum of 1
distarr = reshape(distarr,[length(distarr),1]);
sim = ones(length(distarr),1) ./ distarr;
sim_n = sim / min(sim);
sim_n = (5000/max(sim_n))*sim_n;
sim_r = ceil(sim_n);

weighted_die = [];
for i = 1:length(sim_r)
    weighted_die = [weighted_die i*ones(1,sim_r(i))];
end
weighted_die = weighted_die.';

end









