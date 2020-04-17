%% Example .mat files for loading data
% load mrmrplus_540_150_outstruct.mat
% clearvars -except outstruct
load('mouse_ABA_426LR_brainography.mat');
load('big_ABA_mouse_mats.mat');
load('pervoxLayerMaps.mat');
load('classkey.mat');
load('RajLab_ABAInds.mat');
load('residuals_reg.mat');
% load human_atlas_dk;
% load testgroups.mat;
% load testdata.mat;
% regtotals = outstruct.Bsums;
% normtotals = zeros(426,25);
% regentropy = zeros(1,426);
% entropycalc = @(x) -dot(log(x+eps),x);
% for i = 1:426
%    normtotals(i,:) = regtotals(i,:)/sum(regtotals(i,:)); 
%    regentropy(i) = entropycalc(normtotals(i,:));
% end
% regentropy(isnan(regentropy)) = 0;

%% Field settings for the input structure
%% Anatomical settings
input_struct.voxUreg = 1; 
%Flag setting, 0 for voxel level resolution, 1 for voxels grouped into regions

if input_struct.voxUreg
    q = Q;
    newatlas = zeros(size(renderStruct(1).brain_at));
    for i = 1:length(q)
        f = find(renderStruct(1).brain_at==q(i));
        newatlas(f(:)) = i;
        f2 = find(renderStruct(1).brain_at==(q(i)+213));
        newatlas(f2(:)) = i + 213;
        clear f f2
    end
    renderStruct(1).brain_at = newatlas;
%     datinput = [zeros(30,1);rand(20,1);zeros(20,1);rand(10,1);zeros(6,1)] * 5;
%     datinput = [datinput;datinput]*5;
%     datinput = ones(426,1);
%     datinput = ones(426,1)./regentropy.';
%     datinput(12) = 0; datinput(225) = 0;
    datinput = residuals_reg/mean(residuals_reg);
%     colmap = lines(5);
else    
    testvals = outstruct(3).Bvals(:,18);
    testvals(1:end/2) = 0;
    datamap = zeros(size(GENGDmod));
    for i = 1:length(structList)
        [~,voxinds] = ismember(structIndex{i},nonzerovox);
        allvox = nonzerovox(voxinds);
        curvox = testvals(voxinds);
        datamap(allvox) = curvox;
        clear allvox curvox voxinds
    end
    interpmap = datamap;
    interpmap = imresize3(interpmap,[133 81 115]);
    interpmap(interpmap<0) = 0;
    datinput = interpmap;
    colmap = [flipud((0.2:0.2:1)') ones(5,1) zeros(5,1)];
end

input_struct.brain_atlas = renderStruct(1).brain_at; 
%Setting brain atlas, can be arbitrary species

input_struct.nreg = 426; 
%Setting number of regions

% input_struct.reglabs = listB(:,1); 
%Setting region labels

% testgroups = ones(86,1);
% [sortdat,sortinds] = sort(testdata);
% [~,indresort] = sort(sortinds);
% newgroups = (1:86).';
% newgroups = newgroups(indresort);
% compvec = [newgroups testdata];
% testgroups = newgroups;
% testgroups(69:86) = 2;
% pctiles = prctile(testdata,[25 50 75]);
% testgroups(testdata>pctiles(1) & testdata<=pctiles(2)) = testgroups(testdata>pctiles(1) & testdata<=pctiles(2)) * 2;
% testgroups(testdata>pctiles(2) & testdata<=pctiles(3)) = testgroups(testdata>pctiles(2) & testdata<=pctiles(3)) * 3;
% testgroups(testdata>pctiles(3)) = testgroups(testdata>pctiles(3)) * 4;
% reggroups = randi(5,426,1);
% reggroups = ones(426,1);
% amg = 1:11; cer = 12:22; sub = 23:25; hip = 26:36; hyp = 37:56; neo = 57:94;
% med = 95:119; mid = 120:140; olf = 141:148; pal = 149:156; pns = 157:169;
% str = 170:177; tha = 178:212;
% reggroups = zeros(1,212);
% reggroups(amg) = 1; reggroups(cer) = 2; reggroups(sub) = 3; 
% reggroups(hip) = 4; reggroups(hyp) = 5; reggroups(neo) = 6; 
% reggroups(med) = 7; reggroups(mid) = 8; reggroups(olf) = 9; 
% reggroups(pal) = 10; reggroups(pns) = 11; reggroups(str) = 12; 
% reggroups(tha) = 13;
% reggroups = [reggroups(1:11), 2, reggroups(12:end), reggroups(1:11), 2, reggroups(12:end)];
% input_struct.region_groups = reggroups.'; 
sampleinds = [70,92,186];
newsampleinds = sampleinds+1; newsampleinds = [newsampleinds, newsampleinds+213];
reggroups = ones(426,1)+1; reggroups(newsampleinds) = 1;
input_struct.region_groups = reggroups;
%Example grouping of regions

input_struct.conmat = renderStruct(1).connectivityMatrix;

%% Putting our ISH Data into our Mouse Brain Atlas space
%Operations such as re-indexing and coregsitration should be done outside of the
%function
%Pre-render operations on our data for per voxel rendering

%% Data rendering settings
input_struct.data = datinput*10; 
%Example input data, your data goes here
% input_struct.cmap = flipud(spring(length(unique(testgroups))));
% input_struct.cmap = hsv(length(unique(reggroups)));
input_struct.cmap = [[1 0 0];[0 0 1]];
%Colormap

input_struct.nbin = 1; 
%Number of bins for rendering data based on cumulative distribution, only applicable for per-voxel renderings

input_struct.xfac = 10; 
%Multiplier for number of points rendered per voxel or per region sphere, fractional multipliers not suggested

input_struct.pointsize = 10; 
%Setting size of points rendered

input_struct.size_constant = 500;

input_struct.sphere = 0;

input_struct.centered = [1 2];

%% Image file settings
input_struct.bgcolor = 'w';

input_struct.savenclose = 0; 
%Flag for for 0=opening file in Matlab GUI and 1=saving file and closing .fig

input_struct.img_labels = 'residuals'; 
%Field for setting filename labels, if elected to save

input_struct.img_format = 'tiffn'; 
%Matlab flag for image file format

%% Running the function

brainframe(input_struct)
