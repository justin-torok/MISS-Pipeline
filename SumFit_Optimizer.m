function [sumfitvec,minngvec] = SumFit_Optimizer(fitstruct,method,sigmas,directory)

if nargin < 4
    directory = [cd filesep 'MatFiles']; 
    if nargin < 3
        sigmas = 3000:25:6000;
    end
end

ng_list = zeros(1,length(fitstruct));
ng_param_list = ng_list;
for i = 1:length(fitstruct)
    ng_list(i) = fitstruct(i).nGen;
    if strcmp(method,'MRx3')    
        ng_param_list(i) = fitstruct(i).nGen;
    else
        ng_param_list(i) = fitstruct(i).nG_param;
    end
end
negloglikelihoods = zeros(length(ng_list),length(sigmas));
neglogpriors = negloglikelihoods;
for k = 1:length(sigmas)
    for i = 1:length(fitstruct)
        negloglikelihoods(i,k) = fitstruct(i).negloglikelihood;
        neglogpriors(i,k) = (ng_list(i)^2)/(2*sigmas(k)^2);
    end
end
neglogposteriors = negloglikelihoods + neglogpriors;
mininds = zeros(1,length(sigmas));
sumfitvec = mininds; 
minngvec = mininds; minngparamvec = mininds;
for k = 1:size(mininds,2)
    [~,minind] = min(neglogposteriors(:,k));
    mininds(k) = minind;
    minngvec(k) = ng_list(minind);
    minngparamvec(k) = ng_param_list(minind);
end

load([directory filesep 'PresetInputs.mat'],'regvgene','genevct','gene_names');


outstruct = struct;
if strcmp(method,'MRx3')
    preloadinds = MRx3_Selector(genevct,regvgene,max(minngvec),fitstruct(1).lambda);
else
    preloadinds = [];
end
for k = 1:length(sigmas)
    fprintf('Calculating sum fit, nG parameter value %d/%d\n',k,length(sigmas))

    % Nonnegative matrix inversion
    [E_red,C_red] = GeneSelector(genevct,regvgene,gene_names,minngparamvec(k),...
        fitstruct(1).lambda,method,preloadinds); 
    B = CellDensityInference(E_red,C_red);
    outstruct(k).lambda = fitstruct(1).lambda;
    outstruct(k).nGen = minngvec(k);

    % Density --> Counts correction
    Bcorrected = Density_to_Counts(B,directory);
    outstruct(k).corrB = Bcorrected;

    % Binning voxels into CCF regions
    [sumB,meanB] = Voxel_To_Region(Bcorrected,directory);
    outstruct(k).Bsums = sumB; % total cells per region
    outstruct(k).Bmeans = meanB; % mean cell count per region

    % Lin R Calculation
    LinRstruct = CorrelationsCalc(outstruct,k,directory);
    LinRnames = fieldnames(LinRstruct);
    for n = 1:length(LinRnames)
        curparam_LinR(n,:) = LinRstruct.(LinRnames{n});
    end
    LinR_pv = curparam_LinR(1,end);
    LinR_sst = curparam_LinR(2,end);
    LinR_vip = curparam_LinR(3,end);
    LinR_micro = curparam_LinR(5,end);

    % Calculate adjusted Kendall's tau for layer-type glutamatergic neurons
    ranks = [1 2 3 3 4 4 4];
    cell_inds = 9:15;
    cell_names = {'L2n3','L4','L5IT','L5PT','L6CT','L6IT','L6b'};
    taustruct = TauCalc(outstruct,k,cell_names,cell_inds,ranks,directory);

    % Calculate SumFit criterion
    sumfit = taustruct.tau + LinR_pv + LinR_sst + LinR_vip + LinR_micro;
    sumfitvec(k) = sumfit;
end 

end