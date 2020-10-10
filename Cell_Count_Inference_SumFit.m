function outstruct = Cell_Count_Inference_SumFit(method,ngparamlist,lambda,directory)

if nargin < 4
    directory = [cd filesep 'MatFiles'];
    if nargin < 3
        lambda = 250;
        if nargin < 2
            ngparamlist = [50:50:300,301:600,610:10:1600,1650:100:3750];
            if nargin < 1
                method = 'MRx3';
            end
        end
    end
end

load('Tasic_Inputs.mat','genevct','voxvgene','gene_names');
outstruct = struct;
if strcmp(method,'MRx3')
    preloadinds = MRx3_Selector(genevct,voxvgene,size(voxvgene,2),lambda);
else
    preloadinds = [];
end

for k = 1:length(ngparamlist)
    fprintf('Calculating sum fit, nG parameter value %d/%d\n',k,length(ngparamlist))

    % Nonnegative matrix inversion
    [E_red,C_red] = GeneSelector(genevct,voxvgene,gene_names,ngparamlist(k),...
       lambda,method,preloadinds); 
    B = CellDensityInference(E_red,C_red);
    outstruct(k).lambda = lambda;
    outstruct(k).nGen = ngparamlist(k);

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
    LinR_pv = curparam_LinR(1,end);
    LinR_sst = curparam_LinR(2,end);
    LinR_vip = curparam_LinR(3,end);
    LinR_micro = curparam_LinR(5,end);

    % Calculate adjusted Kendall's tau for layer-type glutamatergic neurons
    ranks = [1 2 3 3 4 4 4];
    cell_inds = 9:15;
    cell_names = {'L2n3','L4','L5IT','L5PT','L6CT','L6IT','L6b'};
    taustruct = TauCalc(outstruct,k,cell_names,cell_inds,ranks,directory);
    outstruct(k).tau = taustruct.tau;

    % Calculate SumFit criterion
    sumfit = taustruct.tau + LinR_pv + LinR_sst + LinR_vip + LinR_micro;
%     sumfitvec(k) = sumfit;
    outstruct(k).sumfit = sumfit;
end 
end
