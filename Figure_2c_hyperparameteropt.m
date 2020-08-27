function [sumfitmat,outstruct] = Figure_2c_hyperparameteropt(lambdas,sigmas,savenclose,directory)

if nargin < 4
    directory = [cd filesep 'MatFiles']; 
    if nargin < 3
        savenclose = 0;
        if nargin < 2
            sigmas = [4000,4100,4200,4300,4400,4500,4600,4700,4800];
            if nargin < 1
                lambdas = 0:50:350; % lambda can only take these values
            end
        end
    end
end

ngcell = cell(1,length(lambdas));
mats2load = ngcell;
for j = 1:length(lambdas)
    mats2load{j} = sprintf('lambda%d_gmm_cv5.mat',lambdas(j));
end

for k = 1:length(sigmas)
    for j = 1:length(mats2load)
        load([directory filesep mats2load{j}],'outstruct','ng_param_list');
        ngcell{j} = ng_param_list;
        for i = 1:length(outstruct)
            if outstruct(i).lambda == 150
                negloglikelihoods(i,j,k) = -(outstruct(i).posteriors^2);
            else
                negloglikelihoods(i,j,k) = -(outstruct(i).likelihood^2);
            end
            neglogpriors(i,j,k) = (ng_param_list(i)^2)/(2*sigmas(k)^2);
        end
        clear outstruct
    end
end

neglogposteriors = negloglikelihoods + neglogpriors;
mininds = zeros(length(mats2load),length(sigmas));
sumfitmat = mininds; 
minngmat = mininds;
for j = 1:size(mininds,1)
    for k = 1:size(mininds,2)
        [~,minind] = min(neglogposteriors(:,j,k));
        mininds(j,k) = minind;
        minngmat(j,k) = ngcell{j}(minind);
    end
end

load([directory filesep 'PresetInputs.mat'],'regvgene','genevct','gene_names');

for j = 1:length(lambdas) 
    outstruct = struct;
    preloadinds = MRx3_Selector(genevct,regvgene,max(minngmat(j,:)),lambdas(j));
    for k = 1:length(sigmas)
        fprintf('Calculating sum fit, nG parameter value %d/%d\n',...
            (k+(j-1)*length(sigmas)),(length(sigmas)*length(lambdas)))
        
        % Nonnegative matrix inversion
        [E_red,C_red] = GeneSelector(genevct,regvgene,gene_names,minngmat(j,k),...
            lambdas(j),'MRx3',preloadinds); 
        B = CellDensityInference(E_red,C_red);
        outstruct(k).lambda = lambdas(j);
        outstruct(k).nGen = minngmat(j,k);
        
        % Density --> Counts correction
        Bcorrected = Density_to_Counts(B,directory);
        outstruct(k).corrB = Bcorrected;
        
        % Binning voxels into CCF regions
        [sumB,meanB] = Voxel_To_Region(Bcorrected,directory);
        outstruct(k).Bsums = sumB; % total cells per region
        outstruct(k).Bmeans = meanB; % mean cell count per region
        
        % Lin R Calculation
        LinRstruct = CorrelationsCalc(outstruct,k,directory);
        outstruct(k).LinR = LinRstruct;
        LinRnames = fieldnames(LinRstruct);
        for n = 1:length(LinRnames)
            curparam_LinR(n,:) = LinRstruct.(LinRnames{n});
        end
        LinR_pv = curparam_LinR(1,3);
        LinR_sst = curparam_LinR(2,3);
        LinR_vip = curparam_LinR(3,3);
        LinR_micro = curparam_LinR(5,3);

        % Calculate adjusted Kendall's tau for layer-type glutamatergic neurons
        ranks = [1 2 3 3 4 4 4];
        cell_inds = 9:15;
        cell_names = {'L2n3','L4','L5IT','L5PT','L6CT','L6IT','L6b'};
        taustruct = TauCalc(outstruct,k,cell_names,cell_inds,ranks,directory);

        % Calculate SumFit criterion
        sumfit = taustruct.tau + LinR_pv + LinR_sst + LinR_vip + LinR_micro;
        sumfitmat(j,k) = sumfit;
    end 
end

% Surface Plotting
% [x,y] = meshgrid(lambdas,sigmas);
% figure; hold on;
% s = surface(x,y,sumfitmat,'FaceAlpha',0.5); colormap hsv;
% s.EdgeColor = 'none';
% xlabel('\lambda'); ylabel('\sigma'); zlabel('\Sigma_{fit}');
% title('Hyperparameter Optimization');
% set(gca,'FontSize',20);

if savenclose
    print('hyperparamopt','-dtiffn');
    close
end
end
