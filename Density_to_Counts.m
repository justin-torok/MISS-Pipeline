function [B_corrected,Bfactor] = Density_to_Counts(B,directory)

if nargin < 2
    directory = [cd filesep 'MatFiles'];
end
neoinds = 57:94;
load([directory filesep 'regionlabs.mat'],'regionlabs');
isneo = ismember(regionlabs,neoinds);
murakami_neodens = 127870; % Total cell density from Murakami et al 2018, cells/mm^3
volfactor = (1/0.2)^3; % Ratio of Murakami to Lein voxel volumes
meanBneo = mean(sum(B(isneo,:),2));
Bfactor = murakami_neodens / (meanBneo * volfactor);
B_corrected = B * Bfactor;


end