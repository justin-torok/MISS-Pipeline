function [B_corrected,Bfactor] = Density_to_Counts(B)

neoinds = 57:94;
load('regionlabs.mat');
isneo = ismember(regionlabs,neoinds);
murakami_neodens = 127870; % Total cell density from Murakami et al 2018, cells/mm^3
% murakami_neodens = 78000000/50246;
volfactor = (1/0.2)^3; % Ratio of Murakami to Lein voxel volumes
% volfactor = 1;
meanBneo = mean(sum(B(isneo,:),2));
% meanBneo = mean(sum(B,2));
Bfactor = murakami_neodens / (meanBneo * volfactor);
B_corrected = B * Bfactor;


end