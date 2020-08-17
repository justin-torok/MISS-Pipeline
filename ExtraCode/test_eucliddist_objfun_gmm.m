%do stuff from MISS_demo here
matdir = '/Users/justintorok/Documents/MATLAB/CellType/MatFiles/MISS';

load([matdir filesep 'PresetInputs.mat'],'regvgene','genevct','classkey','gene_names','ct_labvec','C_indivcells');
% clearvars -except regvgene genevct classkey gene_names ct_labvec C_indivcells matdir
method = 'MRx3';
testnG = [50:1200 1255:100:3755];
% testnG = 3755;
% testnG = [610:10:1600 1650:50:2500 2600:100:3700 3849];
% testnG = [250 350 450 550];
% testnG = [50:50:450 460:10:1600 1650:50:2500 2600:100:3700 3855];
costfun = 'GMMposteriors'; %options: SVM, recluster, GMMposteriors
k = 5;
ng_param_list = testnG;
voxvgene = regvgene;
lambda = 300;
savegroups = 0;
CV = 1;
%     testlambda = 0:50:500;

sumfit_vec = zeros(1,length(ng_param_list));
fprintf('Initializing preloaded gene indices\n')
if lambda == 150 && strcmp(method,'MRx3')
    tempstruct = load([matdir filesep 'MRx3_l150_geneinds.mat']);
    preloadinds = tempstruct.mrmrinds;
    clear tempstruct
elseif (length(ng_param_list) > 5 || strcmp(costfun,'SVM')) && strcmp(method,'MRx3') % heuristic criterion
    preloadinds = MRx3_Selector(genevct,voxvgene,length(gene_names),lambda);
elseif (length(ng_param_list) > 5 || strcmp(costfun,'recluster')) && strcmp(method,'MRx3')
    preloadinds = MRx3_Selector(genevct,voxvgene,length(gene_names),lambda);
else
    preloadinds = [];
end
% preloadinds = mrmrinds;
outstruct = struct;
for i = 1:length(ng_param_list)
    tic
    fprintf('Determining subset, nG parameter value %d/%d\n',i,length(ng_param_list))
    param = ng_param_list(i);
    % Create reduced versions of voxvgene and genevct according to the
    % method and parameter specified by the user
    if strcmp(costfun,'SumFit')
        [E_red,C_red,nGen] = GeneSelector(genevct,voxvgene,gene_names,param,lambda,method,preloadinds);
    else
        [E_red,C_red,nGen,~,C_ind_red] = GeneSelector_IndivCells(genevct,voxvgene,C_indivcells,gene_names,param,lambda,method,preloadinds);
    end
    
    % Infer cell density per voxel in arbitrary units
%     fprintf('Nonnegative matrix inversion, nG parameter value %d/%d\n',i,length(ng_param_list))
%     B = CellDensityInference(E_red,C_red);
%     outstruct(i).Bvals = B; 
    outstruct(i).nGen = nGen;
    if strcmp(method,'MRx3')
        outstruct(i).lambda = lambda;
    end

    % Convert arbitrary densities to counts per voxel
%     Bcorrected = Density_to_Counts(B,matdir);
%     outstruct(i).corrB = Bcorrected;

    % Sum and average over CCF regions
%     [sumB,meanB] = Voxel_To_Region(Bcorrected,matdir);
%     outstruct(i).Bsums = sumB; % total cells per region
%     outstruct(i).Bmeans = meanB; % mean cell count per region

    % Calculate Pearson and Lin R
%     [LinRstruct,PearsonStruct] = CorrelationsCalc(outstruct,i,matdir);
%     outstruct(i).LinR = LinRstruct;
%     outstruct(i).Pearson = PearsonStruct;
%     LinRnames = fieldnames(LinRstruct);
%     Pnames = fieldnames(PearsonStruct);
%     for j = 1:length(LinRnames)
%         curparam_LinR(j,:) = LinRstruct.(LinRnames{j});
%         curparam_Rval(j,:) = PearsonStruct.(Pnames{j});
%     end
%     LinR_pv(i,:) = curparam_LinR(1,:);
%     LinR_sst(i,:) = curparam_LinR(2,:);
%     LinR_vip(i,:) = curparam_LinR(3,:);
%     LinR_all(i,:) = curparam_LinR(4,:);
%     LinR_micro(i,:) = curparam_LinR(5,:);
%     LinR_neuron(i,:) = curparam_LinR(6,:);
%     Rval_pv(i,:) = curparam_Rval(1,:);
%     Rval_sst(i,:) = curparam_Rval(2,:);
%     Rval_vip(i,:) = curparam_Rval(3,:);
%     Rval_all(i,:) = curparam_Rval(4,:);
%     Rval_micro(i,:) = curparam_Rval(5,:);
%     Rval_neuron(i,:) = curparam_Rval(6,:);

    % Calculate adjusted Kendall's tau for layer-type glutamatergic neurons
%     fprintf('Calculating tau, nG parameter value %d/%d\n',i,length(ng_param_list))
%     ranks = [1 2 3 3 4 4 4];
%     cell_inds = 9:15;
%     cell_names = {'L2n3','L4','L5IT','L5PT','L6CT','L6IT','L6b'};
%     taustruct = TauCalc(outstruct,i,cell_names,cell_inds,ranks,matdir);

    % Calculate SumFit criterion
%     fprintf('Calculating sum fit, nG parameter value %d/%d\n',i,length(ng_param_list))
%     sumfit = taustruct.tau + LinR_pv(i,3) + LinR_sst(i,3) + LinR_vip(i,3) + LinR_micro(i,3);
%     sumfit_vec(i) = sumfit;
%     outstruct(i).sumfit = sumfit;
    
    % Calculate classification error (SVM only)
    if strcmp(costfun,'SVM')
        fprintf('Determining SVM classification error, nG parameter value %d/%d\n',i,length(ng_param_list))
        x = (C_ind_red - mean(C_ind_red, 1) ./ std(C_ind_red, 1))';
        y = ct_labvec.';
        options = statset('UseParallel',true);
        Mdl = fitcecoc(x,y,'Options',options);
        rng(0); % set seed
        CVMdl = crossval(Mdl, 'Options', options, 'Kfold', k);
        error_vec(i) = kfoldLoss(CVMdl, 'Options', options);
        outstruct(i).error = error_vec(i);
    elseif strcmp(costfun,'recluster')
        for r = 0:k-1
            rng(r);
            centroids = zeros(size(C_ind_red,1),length(classkey));
            C_testcells = [];
            test_labvec = [];
            for n = 1:length(classkey)
                curcellinds = find(ct_labvec==n);
                testinds = randi(length(curcellinds),1,ceil(length(curcellinds)/k));
                traininds = curcellinds;
                traininds(testinds) = [];
                centroids(:,n) = mean(C_ind_red(:,traininds),2);
                C_testcells = [C_testcells C_ind_red(:,curcellinds(testinds))];
                test_labvec = [test_labvec ct_labvec(curcellinds(testinds))];
            end
            C_testcells = C_testcells.';
            test_labvec = test_labvec.';
            centroids = centroids.';
            C_test = [centroids;C_testcells];
            alldists = pdist(C_test,'euclidean');
            alldists = squareform(alldists);
            test_dists = alldists(26:end,1:25);
            [~,sortinds] = sort(test_dists,2);
            groups_ids = sortinds(:,1);
            if savegroups
                group_mat(:,r+1) = groups_ids;
                test_mat(:,r+1) = test_labvec;
            else
                group_mat = [];
                test_mat = [];
            end
            accuracy(r+1) = sum(test_labvec==groups_ids) / length(test_labvec);
            error(r+1) = 1 - accuracy(r+1);
        end
        outstruct(i).std = std(error);
        outstruct(i).error = mean(error);
        outstruct(i).accuracy = mean(accuracy);
        outstruct(i).testgroups = group_mat;
        outstruct(i).truthgroups = test_mat;
    elseif strcmp(costfun,'GMMposteriors') && CV == 1
        fprintf('Fitting GMM, nG parameter value %d/%d\n',i,length(ng_param_list))
        for r = 0:k-1
            rng(r);
            centroids = zeros(size(C_ind_red,1),length(classkey));
            C_testcells = [];
            test_labvec = [];
            C_traincells = [];
            for n = 1:length(classkey)
                curcellinds = find(ct_labvec==n);
                testinds = randi(length(curcellinds),1,ceil(length(curcellinds)/k));
                traininds = curcellinds;
                traininds(testinds) = [];
                centroids(:,n) = mean(C_ind_red(:,traininds),2);
%                 covmat(:,:,n) = cov(C_ind_red(:,traininds).');
                C_traincells = [C_traincells C_ind_red(:,traininds)];
                C_testcells = [C_testcells C_ind_red(:,curcellinds(testinds))];
                test_labvec = [test_labvec ct_labvec(curcellinds(testinds))];
            end
            C_testcells = C_testcells.';
            C_traincells = C_traincells.';
            covmat = cov(C_traincells);
            test_labvec = test_labvec.';
            centroids = centroids.';
            gm = gmdistribution(centroids,covmat);
%             covmat = gm.Sigma;
            posteriors = posterior(gm,C_testcells);
            test_dists = posteriors;
            [~,sortinds] = sort(test_dists,2,'descend');
            groups_ids = sortinds(:,1);
            indinds = sub2ind(size(posteriors),(1:size(posteriors,1)).',test_labvec);
%             post_probs(:,r+1) = posteriors(indinds);
            post_probs = posteriors(indinds);
            if savegroups
                group_mat(:,r+1) = groups_ids;
                test_mat(:,r+1) = test_labvec;
            else
                group_mat = [];
                test_mat = [];
            end
            accuracy(r+1) = sum(test_labvec==groups_ids) / length(test_labvec);
            error(r+1) = 1 - accuracy(r+1);
            postprob(r+1) = mean(post_probs);
        end
%         postprob = post_probs(:);
        outstruct(i).std_err = std(error);
        outstruct(i).error = mean(error);
        outstruct(i).accuracy = mean(accuracy);
        outstruct(i).likelihood = mean(postprob);
        outstruct(i).std_likelihood = std(postprob);
        outstruct(i).testgroups = group_mat;
        outstruct(i).truthgroups = test_mat;
    elseif strcmp(costfun,'GMMposteriors') && CV == 0
        fprintf('Fitting GMM, nG parameter value %d/%d\n',i,length(ng_param_list))
        rng(0);
        centroids = zeros(size(C_ind_red,1),length(classkey));
        C_testcells = [];
        test_labvec = [];
        C_traincells = [];
        for n = 1:length(classkey)
            curcellinds = find(ct_labvec==n);
            testinds = randi(length(curcellinds),1,ceil(length(curcellinds)/k));
            traininds = curcellinds;
            traininds(testinds) = [];
            centroids(:,n) = mean(C_ind_red(:,traininds),2);
            C_traincells = [C_traincells C_ind_red(:,traininds)];
            C_testcells = [C_testcells C_ind_red(:,curcellinds(testinds))];
            test_labvec = [test_labvec ct_labvec(curcellinds(testinds))];
        end
        C_testcells = C_testcells.';
        C_traincells = C_traincells.';
        covmat = cov(C_traincells);
        test_labvec = test_labvec.';
        centroids = centroids.';
        gm = gmdistribution(centroids,covmat);
        posteriors = posterior(gm,C_testcells);
        test_dists = posteriors;
        [~,sortinds] = sort(test_dists,2,'descend');
        groups_ids = sortinds(:,1);
        indinds = sub2ind(size(posteriors),(1:size(posteriors,1)).',test_labvec);
        post_probs = posteriors(indinds);
        if savegroups
            group_mat = groups_ids;
            test_mat = test_labvec;
        else
            group_mat = [];
            test_mat = [];
        end
        accuracy = sum(test_labvec==groups_ids) / length(test_labvec);
        error = 1 - accuracy;
        outstruct(i).std_err = std(error);
        outstruct(i).error = error;
        outstruct(i).accuracy = accuracy;
        outstruct(i).posteriors = post_probs;
        outstruct(i).std_post = std(1-post_probs);
        outstruct(i).testgroups = group_mat;
        outstruct(i).truthgroups = test_mat;
    end
    fprintf('Done, nG parameter value %d/%d\n',i,length(ng_param_list))
    toc
end

% outstruct = outstruct(1:end-1);
% ng_param_list = ng_param_list(1:end-1);
% clear posteriors
% for i = 1:length(outstruct)
%     error(i) = outstruct(i).error;
%     accuracy(i) = outstruct(i).accuracy;
%     posteriors(i) = outstruct(i).posteriors;
% end
% norm_error = (error - min(error)) / (max(error) - min(error));
% post_notright = 1 - posteriors;
% norm_post = (post_notright - min(post_notright)) / (max(post_notright) - min(post_notright));
% normnG = ((ng_param_list - min(ng_param_list)) / (3855 - min(ng_param_list))).';
% normnG = normnG.';
% dist2origin = sqrt((normnG-0).^2 + (norm_error-0).^2);
% [~,elbowind] = min(dist2origin);
% 
% figure;
% plot(normnG,norm_error,'Color',[0.9 0.4 0.2],'LineWidth',3); hold on;
% plot([normnG(elbowind) normnG(elbowind)],[0 1],'Color',[0 0 0],'LineWidth',2,'LineStyle','--'); hold on;
% title('Recluster Error vs. nG Elbow Curve','FontSize',18);
% xlabel('Norm. nG','FontSize',18);
% ylabel('Norm. Error','FontSize',18);
% xticks([0 0.25 0.5 0.75 1]);
% yticks([0 0.25 0.5 0.75 1]);
% set(gca,'FontSize',16);
% 
% figure;
% plot(normnG,norm_post,'Color',[0.9 0.4 0.2],'LineWidth',3); hold on;
% plot([normnG(elbowind) normnG(elbowind)],[0 1],'Color',[0 0 0],'LineWidth',2,'LineStyle','--'); hold on;
% title('1-Posteriors vs. nG Elbow Curve','FontSize',18);
% xlabel('Norm. nG','FontSize',18);
% ylabel('Norm. 1-Posteriors','FontSize',18);
% xticks([0 0.25 0.5 0.75 1]);
% yticks([0 0.25 0.5 0.75 1]);
% set(gca,'FontSize',16);
% 
save([matdir filesep 'lambda' num2str(lambda) '_gmm_cv5.mat'],...
'accuracy','classkey','costfun','error','lambda','matdir',...
'method','ng_param_list','outstruct','preloadinds')