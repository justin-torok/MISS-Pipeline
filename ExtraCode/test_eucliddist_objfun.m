%do stuff from MISS_demo here
matdir = '/Users/justintorok/Documents/MATLAB/CellType/MatFiles/MISS';

load([matdir filesep 'Tasic_Inputs.mat'],'voxvgene','genevct','classkey','gene_names','ct_labvec','C_indivcells');
% clearvars -except regvgene genevct classkey gene_names ct_labvec C_indivcells matdir
% load([matdir filesep 'MRx3_l150_geneinds.mat']);
method = 'MRx3';
% testnG = [50:20:300 301:600 620:20:800 800:50:1200 1255:100:3755];
testnG = 300:20:800;
% testnG = 50;
% testnG = [610:10:1600 1650:50:2500 2600:100:3700 3849];
% testnG = [250 350 450 550];
% testnG = [50:50:450 460:10:1600 1650:50:2500 2600:100:3700 3855];
costfun = 'SumFit'; %options: SVM, recluster, GMMposteriors
resflag_ = '';
% k = 5;
ng_param_list = testnG;
lambda = 250;
savegroups = 0;
% CV = 0;
% dense_factor = '';
fname_label = 'l250_ng529_lsqlin';
%     testlambda = 0:50:500;

% preloadinds = MRx3_Selector_TrueF(genevct,voxvgene,length(gene_names),lambda,C_indivcells,ct_labvec);
% sumfit_vec = zeros(1,length(ng_param_list));
fprintf('Initializing preloaded gene indices\n')
if (length(ng_param_list) > 5 || strcmp(costfun,'SVM')) && strcmp(method,'MRx3') % heuristic criterion
    preloadinds = MRx3_Selector(genevct,voxvgene,length(gene_names),lambda);
elseif (length(ng_param_list) > 5 || strcmp(costfun,'recluster')) && strcmp(method,'MRx3')
    preloadinds = MRx3_Selector(genevct,voxvgene,length(gene_names),lambda);
else
    preloadinds = [];
end
% preloadinds = mrmrinds;
% preloadinds = [];

tauvec = zeros(length(ng_param_list),1);
LinR_pv = zeros(length(ng_param_list),4);
LinR_sst = zeros(length(ng_param_list),4);
LinR_vip = zeros(length(ng_param_list),4);
LinR_micro = zeros(length(ng_param_list),4);
sumfit_vec = zeros(length(ng_param_list),1);
for i = 1:length(ng_param_list)
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
    fprintf('Nonnegative matrix inversion, nG parameter value %d/%d\n',i,length(ng_param_list))
    [B,rvals] = CellDensityInference_Emod(E_red,C_red,resflag_);
%     B = B.* rvals;
%     B = CellDensityInference(E_red,C_red);
    outstruct(i).Bvals = B; 
    outstruct(i).nGen = nGen;
    if strcmp(method,'MRx3')
        outstruct(i).lambda = lambda;
    end

    % Convert arbitrary densities to counts per voxel
    Bcorrected = Density_to_Counts(B,matdir);
    outstruct(i).corrB = Bcorrected;

    % Sum and average over CCF regions
    [sumB,meanB] = Voxel_To_Region(Bcorrected,matdir);
    outstruct(i).Bsums = sumB; % total cells per region
    outstruct(i).Bmeans = meanB; % mean cell count per region

    % Calculate Pearson and Lin R
    [LinRstruct,PearsonStruct] = CorrelationsCalc(outstruct,i,matdir);
    outstruct(i).LinR = LinRstruct;
    outstruct(i).Pearson = PearsonStruct;
    LinRnames = fieldnames(LinRstruct);
%     Pnames = fieldnames(PearsonStruct);
    for j = 1:length(LinRnames)
        curparam_LinR(j,:) = LinRstruct.(LinRnames{j});
%         curparam_Rval(j,:) = PearsonStruct.(Pnames{j});
    end
    LinR_pv(i,:) = curparam_LinR(1,:);
    LinR_sst(i,:) = curparam_LinR(2,:);
    LinR_vip(i,:) = curparam_LinR(3,:);
%     LinR_all(i,:) = curparam_LinR(4,:);
    LinR_micro(i,:) = curparam_LinR(5,:);
%     LinR_neuron(i,:) = curparam_LinR(6,:);
%     Rval_pv(i,:) = curparam_Rval(1,:);
%     Rval_sst(i,:) = curparam_Rval(2,:);
%     Rval_vip(i,:) = curparam_Rval(3,:);
%     Rval_all(i,:) = curparam_Rval(4,:);
%     Rval_micro(i,:) = curparam_Rval(5,:);
%     Rval_neuron(i,:) = curparam_Rval(6,:);

    % Calculate adjusted Kendall's tau for layer-type glutamatergic neurons
    fprintf('Calculating tau, nG parameter value %d/%d\n',i,length(ng_param_list))
    ranks = [1 2 3 3 4 4 4];
    cell_inds = 9:15;
    cell_names = {'L2n3','L4','L5IT','L5PT','L6CT','L6IT','L6b'};
    taustruct = TauCalc(outstruct,i,cell_names,cell_inds,ranks,matdir);
    tauvec(i) = taustruct.tau;
    outstruct(i).tau = taustruct.tau;

    % Calculate SumFit criterion
    fprintf('Calculating sum fit, nG parameter value %d/%d\n',i,length(ng_param_list))
    sumfit = taustruct.tau + LinR_pv(i,end) + LinR_sst(i,end) + LinR_vip(i,end) + LinR_micro(i,end);
    sumfit_vec(i) = sumfit;
    outstruct(i).sumfit = sumfit;
    
    % Calculate classification error (SVM only)
%     if strcmp(costfun,'SVM')
%         fprintf('Determining SVM classification error, nG parameter value %d/%d\n',i,length(ng_param_list))
%         x = (C_ind_red - mean(C_ind_red, 1) ./ std(C_ind_red, 1))';
%         y = ct_labvec.';
%         options = statset('UseParallel',true);
%         Mdl = fitcecoc(x,y,'Options',options);
%         rng(0); % set seed
%         CVMdl = crossval(Mdl, 'Options', options, 'Kfold', k);
%         error_vec(i) = kfoldLoss(CVMdl, 'Options', options);
%         outstruct(i).error = error_vec(i);
%     elseif strcmp(costfun,'recluster')
%         for r = 0:k-1
%             rng(r);
%             centroids = zeros(size(C_ind_red,1),length(classkey));
%             C_testcells = [];
%             test_labvec = [];
%             for n = 1:length(classkey)
%                 curcellinds = find(ct_labvec==n);
%                 testinds = randi(length(curcellinds),1,ceil(length(curcellinds)/k));
%                 traininds = curcellinds;
%                 traininds(testinds) = [];
%                 centroids(:,n) = mean(C_ind_red(:,traininds),2);
%                 C_testcells = [C_testcells C_ind_red(:,curcellinds(testinds))];
%                 test_labvec = [test_labvec ct_labvec(curcellinds(testinds))];
%             end
%             C_testcells = C_testcells.';
%             test_labvec = test_labvec.';
%             centroids = centroids.';
%             C_test = [centroids;C_testcells];
%             alldists = pdist(C_test,'euclidean');
%             alldists = squareform(alldists);
%             test_dists = alldists(26:end,1:25);
%             [~,sortinds] = sort(test_dists,2);
%             groups_ids = sortinds(:,1);
%             if savegroups
%                 group_mat(:,r+1) = groups_ids;
%                 test_mat(:,r+1) = test_labvec;
%             else
%                 group_mat = [];
%                 test_mat = [];
%             end
%             accuracy(r+1) = sum(test_labvec==groups_ids) / length(test_labvec);
%             error(r+1) = 1 - accuracy(r+1);
%         end
%         outstruct(i).std = std(error);
%         outstruct(i).error = mean(error);
%         outstruct(i).accuracy = mean(accuracy);
%         outstruct(i).testgroups = group_mat;
%         outstruct(i).truthgroups = test_mat;
%     elseif strcmp(costfun,'GMMposteriors') && CV == 1
%         for r = 0:k-1
%             rng(r);
%             centroids = zeros(size(C_ind_red,1),length(classkey));
%             C_testcells = [];
%             test_labvec = [];
%             C_traincells = [];
%             for n = 1:length(classkey)
%                 curcellinds = find(ct_labvec==n);
%                 testinds = randi(length(curcellinds),1,ceil(length(curcellinds)/k));
%                 traininds = curcellinds;
%                 traininds(testinds) = [];
%                 centroids(:,n) = mean(C_ind_red(:,traininds),2);
% %                 covmat(:,:,n) = cov(C_ind_red(:,traininds).');
%                 C_traincells = [C_traincells C_ind_red(:,traininds)];
%                 C_testcells = [C_testcells C_ind_red(:,curcellinds(testinds))];
%                 test_labvec = [test_labvec ct_labvec(curcellinds(testinds))];
%             end
%             C_testcells = C_testcells.';
%             C_traincells = C_traincells.';
%             covmat = cov(C_traincells);
%             test_labvec = test_labvec.';
%             centroids = centroids.';
%             gm = gmdistribution(centroids,covmat);
% %             covmat = gm.Sigma;
%             posteriors = posterior(gm,C_testcells);
%             test_dists = posteriors;
%             [~,sortinds] = sort(test_dists,2,'descend');
%             groups_ids = sortinds(:,1);
%             indinds = sub2ind(size(posteriors),(1:size(posteriors,1)).',test_labvec);
%             post_probs(:,r+1) = posteriors(indinds);
%             if savegroups
%                 group_mat(:,r+1) = groups_ids;
%                 test_mat(:,r+1) = test_labvec;
%             else
%                 group_mat = [];
%                 test_mat = [];
%             end
%             accuracy(r+1) = sum(test_labvec==groups_ids) / length(test_labvec);
%             error(r+1) = 1 - accuracy(r+1);
% %             postprob(r+1) = mean(post_probs);
%         end
%         postprob = post_probs(:);
%         outstruct(i).std_err = std(error);
%         outstruct(i).error = mean(error);
%         outstruct(i).accuracy = mean(accuracy);
%         outstruct(i).posteriors = mean(postprob);
%         outstruct(i).std_post = std(postprob);
%         outstruct(i).testgroups = group_mat;
%         outstruct(i).truthgroups = test_mat;
%     elseif strcmp(costfun,'GMMposteriors') && CV == 0
%         rng(0);
%         centroids = zeros(size(C_ind_red,1),length(classkey));
%         C_testcells = [];
%         test_labvec = [];
%         C_traincells = [];
%         for n = 1:length(classkey)
%             curcellinds = find(ct_labvec==n);
%             testinds = randi(length(curcellinds),1,ceil(length(curcellinds)/k));
%             traininds = curcellinds;
%             traininds(testinds) = [];
%             centroids(:,n) = mean(C_ind_red(:,traininds),2);
%             C_traincells = [C_traincells C_ind_red(:,traininds)];
%             C_testcells = [C_testcells C_ind_red(:,curcellinds(testinds))];
%             test_labvec = [test_labvec ct_labvec(curcellinds(testinds))];
%         end
%         C_testcells = C_testcells.';
%         C_traincells = C_traincells.';
%         covmat = cov(C_traincells);
%         test_labvec = test_labvec.';
%         centroids = centroids.';
%         gm = gmdistribution(centroids,covmat);
%         posteriors = posterior(gm,C_testcells);
%         test_dists = posteriors;
%         [~,sortinds] = sort(test_dists,2,'descend');
%         groups_ids = sortinds(:,1);
%         indinds = sub2ind(size(posteriors),(1:size(posteriors,1)).',test_labvec);
%         post_probs = posteriors(indinds);
%         if savegroups
%             group_mat = groups_ids;
%             test_mat = test_labvec;
%         else
%             group_mat = [];
%             test_mat = [];
%         end
%         accuracy = sum(test_labvec==groups_ids) / length(test_labvec);
%         error = 1 - accuracy;
%         outstruct(i).std_err = std(error);
%         outstruct(i).error = error;
%         outstruct(i).accuracy = accuracy;
%         outstruct(i).posteriors = mean(post_probs);
%         outstruct(i).std_post = std(post_probs);
%         outstruct(i).testgroups = group_mat;
%         outstruct(i).truthgroups = test_mat;
%     end
    fprintf('Done, nG parameter value %d/%d\n',i,length(ng_param_list))
end

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
% figure;
% plot(ng_param_list,sumfit_vec,'Color',[0.9 0.4 0.2],'LineWidth',3); hold on;
% % plot([normnG(elbowind) normnG(elbowind)],[0 1],'Color',[0 0 0],'LineWidth',2,'LineStyle','--'); hold on;
% title('n_G Vs. \Sigma_f_i_t','FontSize',18);
% xlabel(['n_G, max @ ' num2str(ng_param_list(sumfit_vec==max(sumfit_vec)))],'FontSize',18);
% ylabel(['\Sigma_f_i_t, max = ' num2str(max(sumfit_vec))],'FontSize',18);
% % xticks([0 0.25 0.5 0.75 1]);
% % yticks([0 0.25 0.5 0.75 1]);
% set(gca,'FontSize',16);
% 
% figure;
% plot(ng_param_list,LinR_pv(:,end),'Color',[0.9 0.4 0.2],'LineWidth',3); hold on;
% % plot([normnG(elbowind) normnG(elbowind)],[0 1],'Color',[0 0 0],'LineWidth',2,'LineStyle','--'); hold on;
% title('n_G Vs. Pvalb+ R_C','FontSize',18);
% xlabel(['n_G, max @ ' num2str(ng_param_list(LinR_pv(:,end)==max(LinR_pv(:,end))))],'FontSize',18);
% ylabel(['Pvalb+ R_C, max = ' num2str(max(LinR_pv(:,end)))],'FontSize',18);
% % xticks([0 0.25 0.5 0.75 1]);
% % yticks([0 0.25 0.5 0.75 1]);
% set(gca,'FontSize',16);
% 
% figure;
% plot(ng_param_list,LinR_sst(:,end),'Color',[0.9 0.4 0.2],'LineWidth',3); hold on;
% % plot([normnG(elbowind) normnG(elbowind)],[0 1],'Color',[0 0 0],'LineWidth',2,'LineStyle','--'); hold on;
% title('n_G Vs. Sst+ R_C','FontSize',18);
% xlabel(['n_G, max @ ' num2str(ng_param_list(LinR_sst(:,end)==max(LinR_sst(:,end))))],'FontSize',18);
% ylabel(['Sst+ R_C, max = ' num2str(max(LinR_sst(:,end)))],'FontSize',18);
% % xticks([0 0.25 0.5 0.75 1]);
% % yticks([0 0.25 0.5 0.75 1]);
% set(gca,'FontSize',16);
% 
% figure;
% plot(ng_param_list,LinR_vip(:,end),'Color',[0.9 0.4 0.2],'LineWidth',3); hold on;
% % plot([normnG(elbowind) normnG(elbowind)],[0 1],'Color',[0 0 0],'LineWidth',2,'LineStyle','--'); hold on;
% title('n_G Vs. Vip+ R_C','FontSize',18);
% xlabel(['n_G, max @ ' num2str(ng_param_list(LinR_vip(:,end)==max(LinR_vip(:,end))))],'FontSize',18);
% ylabel(['Vip+ R_C, max = ' num2str(max(LinR_vip(:,end)))],'FontSize',18);
% % xticks([0 0.25 0.5 0.75 1]);
% % yticks([0 0.25 0.5 0.75 1]);
% set(gca,'FontSize',16);
% 
% figure;
% plot(ng_param_list,LinR_micro(:,end),'Color',[0.9 0.4 0.2],'LineWidth',3); hold on;
% % plot([normnG(elbowind) normnG(elbowind)],[0 1],'Color',[0 0 0],'LineWidth',2,'LineStyle','--'); hold on;
% title('n_G Vs. Micro R_C','FontSize',18);
% xlabel(['n_G, max @ ' num2str(ng_param_list(LinR_micro(:,end)==max(LinR_micro(:,end))))],'FontSize',18);
% ylabel(['Vip+ R_C, max = ' num2str(max(LinR_micro(:,end)))],'FontSize',18);
% % xticks([0 0.25 0.5 0.75 1]);
% % yticks([0 0.25 0.5 0.75 1]);
% set(gca,'FontSize',16);
% 
% figure;
% plot(ng_param_list,tauvec,'Color',[0.9 0.4 0.2],'LineWidth',3); hold on;
% % plot([normnG(elbowind) normnG(elbowind)],[0 1],'Color',[0 0 0],'LineWidth',2,'LineStyle','--'); hold on;
% title('n_G Vs. \tau_{Adj}','FontSize',18);
% xlabel(['n_G, max @ ' num2str(ng_param_list(tauvec==max(tauvec)))],'FontSize',18);
% ylabel(['\tau_{Adj}, max = ' num2str(max(tauvec))],'FontSize',18);
% % xticks([0 0.25 0.5 0.75 1]);
% % yticks([0 0.25 0.5 0.75 1]);
% set(gca,'FontSize',16);
% 
% sumfit_neo = LinR_pv(:,1) + LinR_sst(:,1) + LinR_vip(:,1) + LinR_micro(:,1) + tauvec;
% sumfit_fore = LinR_pv(:,2) + LinR_sst(:,2) + LinR_vip(:,2) + LinR_micro(:,2) + tauvec;
% 
% figure;
% plot(ng_param_list,sumfit_neo,'Color',[0.9 0.4 0.2],'LineWidth',3); hold on;
% % plot([normnG(elbowind) normnG(elbowind)],[0 1],'Color',[0 0 0],'LineWidth',2,'LineStyle','--'); hold on;
% title('n_G Vs. \Sigma_f_i_t, Nctx.','FontSize',18);
% xlabel(['n_G, max @ ' num2str(ng_param_list(sumfit_neo==max(sumfit_neo)))],'FontSize',18);
% ylabel(['\Sigma_f_i_t, max = ' num2str(max(sumfit_neo))],'FontSize',18);
% % xticks([0 0.25 0.5 0.75 1]);
% % yticks([0 0.25 0.5 0.75 1]);
% set(gca,'FontSize',16);
% 
% figure;
% plot(ng_param_list,LinR_pv(:,1),'Color',[0.9 0.4 0.2],'LineWidth',3); hold on;
% % plot([normnG(elbowind) normnG(elbowind)],[0 1],'Color',[0 0 0],'LineWidth',2,'LineStyle','--'); hold on;
% title('n_G Vs. Pvalb+ R_C, Nctx.','FontSize',18);
% xlabel(['n_G, max @ ' num2str(ng_param_list(LinR_pv(:,1)==max(LinR_pv(:,1))))],'FontSize',18);
% ylabel(['Pvalb+ R_C, max = ' num2str(max(LinR_pv(:,1)))],'FontSize',18);
% % xticks([0 0.25 0.5 0.75 1]);
% % yticks([0 0.25 0.5 0.75 1]);
% set(gca,'FontSize',16);
% 
% figure;
% plot(ng_param_list,LinR_sst(:,1),'Color',[0.9 0.4 0.2],'LineWidth',3); hold on;
% % plot([normnG(elbowind) normnG(elbowind)],[0 1],'Color',[0 0 0],'LineWidth',2,'LineStyle','--'); hold on;
% title('n_G Vs. Sst+ R_C, Nctx.','FontSize',18);
% xlabel(['n_G, max @ ' num2str(ng_param_list(LinR_sst(:,1)==max(LinR_sst(:,1))))],'FontSize',18);
% ylabel(['Sst+ R_C, max = ' num2str(max(LinR_sst(:,1)))],'FontSize',18);
% % xticks([0 0.25 0.5 0.75 1]);
% % yticks([0 0.25 0.5 0.75 1]);
% set(gca,'FontSize',16);
% 
% figure;
% plot(ng_param_list,LinR_vip(:,1),'Color',[0.9 0.4 0.2],'LineWidth',3); hold on;
% % plot([normnG(elbowind) normnG(elbowind)],[0 1],'Color',[0 0 0],'LineWidth',2,'LineStyle','--'); hold on;
% title('n_G Vs. Vip+ R_C, Nctx.','FontSize',18);
% xlabel(['n_G, max @ ' num2str(ng_param_list(LinR_vip(:,1)==max(LinR_vip(:,1))))],'FontSize',18);
% ylabel(['Vip+ R_C, max = ' num2str(max(LinR_vip(:,1)))],'FontSize',18);
% % xticks([0 0.25 0.5 0.75 1]);
% % yticks([0 0.25 0.5 0.75 1]);
% set(gca,'FontSize',16);
% 
% figure;
% plot(ng_param_list,sumfit_fore,'Color',[0.9 0.4 0.2],'LineWidth',3); hold on;
% % plot([normnG(elbowind) normnG(elbowind)],[0 1],'Color',[0 0 0],'LineWidth',2,'LineStyle','--'); hold on;
% title('n_G Vs. \Sigma_f_i_t, F.B.','FontSize',18);
% xlabel(['n_G, max @ ' num2str(ng_param_list(sumfit_fore==max(sumfit_fore)))],'FontSize',18);
% ylabel(['\Sigma_f_i_t, max = ' num2str(max(sumfit_fore))],'FontSize',18);
% % xticks([0 0.25 0.5 0.75 1]);
% % yticks([0 0.25 0.5 0.75 1]);
% set(gca,'FontSize',16);
% 
% figure;
% plot(ng_param_list,LinR_pv(:,2),'Color',[0.9 0.4 0.2],'LineWidth',3); hold on;
% % plot([normnG(elbowind) normnG(elbowind)],[0 1],'Color',[0 0 0],'LineWidth',2,'LineStyle','--'); hold on;
% title('n_G Vs. Pvalb+ R_C, F.B.','FontSize',18);
% xlabel(['n_G, max @ ' num2str(ng_param_list(LinR_pv(:,2)==max(LinR_pv(:,1))))],'FontSize',18);
% ylabel(['Pvalb+ R_C, max = ' num2str(max(LinR_pv(:,2)))],'FontSize',18);
% % xticks([0 0.25 0.5 0.75 1]);
% % yticks([0 0.25 0.5 0.75 1]);
% set(gca,'FontSize',16);
% 
% figure;
% plot(ng_param_list,LinR_sst(:,2),'Color',[0.9 0.4 0.2],'LineWidth',3); hold on;
% % plot([normnG(elbowind) normnG(elbowind)],[0 1],'Color',[0 0 0],'LineWidth',2,'LineStyle','--'); hold on;
% title('n_G Vs. Sst+ R_C, F.B.','FontSize',18);
% xlabel(['n_G, max @ ' num2str(ng_param_list(LinR_sst(:,2)==max(LinR_sst(:,2))))],'FontSize',18);
% ylabel(['Sst+ R_C, max = ' num2str(max(LinR_sst(:,2)))],'FontSize',18);
% % xticks([0 0.25 0.5 0.75 1]);
% % yticks([0 0.25 0.5 0.75 1]);
% set(gca,'FontSize',16);
% 
% figure;
% plot(ng_param_list,LinR_vip(:,2),'Color',[0.9 0.4 0.2],'LineWidth',3); hold on;
% % plot([normnG(elbowind) normnG(elbowind)],[0 1],'Color',[0 0 0],'LineWidth',2,'LineStyle','--'); hold on;
% title('n_G Vs. Vip+ R_C, F.B.','FontSize',18);
% xlabel(['n_G, max @ ' num2str(ng_param_list(LinR_vip(:,2)==max(LinR_vip(:,2))))],'FontSize',18);
% ylabel(['Vip+ R_C, max = ' num2str(max(LinR_vip(:,2)))],'FontSize',18);
% % xticks([0 0.25 0.5 0.75 1]);
% % yticks([0 0.25 0.5 0.75 1]);
% set(gca,'FontSize',16);
% Figure_4b_interneuron_scatter(outstruct,1,0,matdir);
% Figure_4ef_glia(outstruct,1,0,matdir);
% Figure_5ab_taulayerslice(outstruct,1,[25,31,36],0,matdir);
% save([fname_label num2str(lambda) '_sumfit.mat'],'classkey','lambda','matdir','method','ng_param_list','outstruct','preloadinds','LinR_pv','LinR_sst','LinR_vip','LinR_micro','tauvec','sumfit_vec','-v7.3')
% save(['lambda' num2str(lambda) '_superfine_recluster.mat'],'accuracy','classkey','costfun','elbowind','error','lambda','matdir','method','mrmrinds','ng_param_list','outstruct','posteriors','preloadinds')