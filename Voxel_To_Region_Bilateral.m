function [cell_types_sum,cell_types_mean] = Voxel_To_Region_Bilateral(D,directory)

if nargin < 2
    directory = [cd filesep 'MatFiles'];
end

load([directory filesep 'regionlabs.mat'],'regionlabs');

if size(D,1) < size(D,2)
    D = D.';
end

cell_types = zeros(length(unique(regionlabs)),size(D,2));
for i = 1:length(unique(regionlabs))
    indi = regionlabs == i;
    regi = D(indi,:);
    cell_types_sum(i,:) = sum(regi);
    cell_types_mean(i,:) = mean(regi);
end
zvec = zeros(1,size(cell_types_sum,2));
cell_types_sum = [cell_types_sum(1:11,:);zvec;cell_types_sum(12:end,:)];
cell_types_mean = [cell_types_mean(1:11,:);zvec;cell_types_mean(12:end,:)];
end
