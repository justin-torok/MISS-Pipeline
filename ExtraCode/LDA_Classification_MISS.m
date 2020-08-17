function [error_vec, mean_post_vec] = LDA_Classification_MISS(C_ind_red_, ct_labvec_, options_, n_groups)
% Function that calculates SVM classification accuracy using fitcecoc() and
% k-fold cross-validation

if nargin < 4
    n_groups = 10;
    if nargin < 3
        options_ = statset('UseParallel',true); 
    end
end
rng(0); % 0 is the default seed, for reproducibility
if options_.UseParallel
    parpool(2); % 4 workers, standard setting that does not break on powerful machines
end

% Normalizing sample expression matrix for classification
C_ind_red_ = (C_ind_red_ - mean(C_ind_red_, 1) ./ std(C_ind_red_, 1))';

% Defining the indices of the blocks (sub matrices) of C_ind_red
blockcellinds = cell(n_groups, length(unique(ct_labvec_)));
for j = 1:length(unique(ct_labvec_)) 
    j_inds = find(ct_labvec_ == j);
    num_samples = length(j_inds);
    random_sort_inds = randperm(num_samples);
    j_inds_rand = j_inds(random_sort_inds);
    startind = 1;
    for i = 1:n_groups
        if i < n_groups
            groupsize = round(num_samples/n_groups);
        else
            groupsize = num_samples - ((n_groups-1)*round(num_samples/n_groups)) ;
        end
        blockcellinds{i,j} = j_inds_rand(startind:(startind+groupsize-1));
        startind = startind + groupsize;
    end
end
clear i j

% Creating the blocks (indices)
testinds = cell(1,n_groups);
traininds = testinds;
% sub_C_train_cell = cell(1, n_groups);
% sub_ct_train_labvec = sub_C_train_cell;
% sub_C_test_cell = sub_C_train_cell;
% sub_ct_test_labvec = sub_ct_train_labvec;
for i = 1:n_groups
    blockinds_ = [];
%     sub_ct_labvec_ = blockinds_;
    for j = 1:length(unique(ct_labvec_))
        blockinds_ = cat(2,blockinds_,blockcellinds{i,j});
%         sub_ct_labvec_ = cat(2,sub_ct_labvec_,j*ones(1,length(blockcellinds{i,j})));
    end
%     sub_C_test_cell{i} = C_ind_red_(blockinds_,:);
    nonblockinds = setdiff(1:size(C_ind_red_,1), blockinds_);
%     sub_C_train_cell{i} = C_ind_red_(nonblockinds,:);
%     sub_ct_test_labvec{i} = ct_labvec_(blockinds_);
%     sub_ct_train_labvec{i} = ct_labvec_(nonblockinds);
    testinds{i} = blockinds_;
    traininds{i} = nonblockinds;
    clear blockinds_ nonblockinds
end
% clear C_ind_red_ % save RAM, large matrix that is redundant

% SVM classification
error_vec = zeros(1,n_groups);
mean_post_vec = error_vec;
parfor i = 1:n_groups
    x_train = C_ind_red_(traininds{i},:);
    x_test = C_ind_red_(testinds{i},:);
    y_train = ct_labvec_(traininds{i});
    y_test = ct_labvec_(testinds{i});
    [testclass,~,posteriors] = classify(x_test,x_train,y_train,'linear');
    testposteriors = posteriors(sub2ind(size(posteriors),1:size(posteriors,1),y_test));
    mean_post_vec(i) = mean(testposteriors);
    error_vec(i) = 1 - (sum(testclass == y_test.')/length(y_test));
%     x = sub_C_train_cell{i};
%     y = (sub_ct_train_labvec{i}).';
%     CVMdl = fitcecoc(x,y,'Learners','svm','Coding','onevsall','CrossVal','on','Verbose',1,'Options',options_);
%     CVMdl = crossval(Mdl, 'Options', options_, 'Kfold', k_);
%     error_vec(i) = kfoldLoss(CVMdl, 'Options', options_);
end  

if options_.UseParallel
    delete(gcp('nocreate'));
end
end

