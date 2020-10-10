% matdir = '/Users/christophermezias/Documents/MISS_General/MatFiles';
load([matdir filesep 'zeisel_l250_ng1168_outputs.mat']);
load([matdir filesep 'Zeisel_Inputs.mat'],'classkey');
idx = 1;
% T_z = Figure_7a_dendrograms_Zeisel(outstruct,idx,'fore',0,matdir);
% load([matdir filesep 'lambda150_rng0_superfine_300to600.mat']);
% Figure5_AllCT_Scatter(outstruct,1,0,matdir)

% visinds = [161 153 173 156 164 27 83 65 92 47];
visinds = [153 27 83 92 117];


% Cell_Type_Brainframe(outstruct,1,visinds,0,[],'zeisel',matdir)
% slicelocs = [35 40 65];
slicelocs = {[35 45 60];[35 40 65];[40 45 65];[35 45 60];[30 50 65]};
% Figure_4ab_taulayerslice(outstruct,245,slicelocs,0,matdir)
% slicelocs = [30 35 40 45 50 55 60 65];
% Density_Slice_Maps(outstruct,1,visinds,slicelocs,1,'zeisel',matdir)

Figure5_ZeiselCompMaps(outstruct,1,visinds,slicelocs,1,matdir);

% Cell_Type_Brainframe(outstruct,1,visinds,0,[],'zeisel',matdir)

% Density_Slice_Maps_test(outstruct,1,visinds,slicelocs,1,'zeisel',matdir)

% Zeisel_MISS_correlationmaps(visinds,'MRx3',1168,250,slicelocs,1,'zeisel',matdir)
% Figure_2d_correlationmaps(visinds,method,ngen_param,lambda,slicelocs,savenclose,directory)
    