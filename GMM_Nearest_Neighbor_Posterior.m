function gmmstruct = GMM_Nearest_Neighbor_Posterior(C_indivcells_, ct_labvec_, k_, crossval_, savegroups_)

if nargin < 5
    savegroups_ = 0;
    if nargin < 4
        crossval_ = 1;
        if nargin < 3
            k_ = 5;
        end
    end
end
C_ind_red = C_indivcells_;
ct_labvec = ct_labvec_;
gmmstruct = struct;
nts = length(unique(ct_labvec));

if crossval_ == 1
    for r = 0:k_-1
        rng(r);
        centroids = zeros(size(C_ind_red,1),nts);
        C_testcells = [];
        test_labvec = [];
        C_traincells = [];
        for n = 1:nts
            curcellinds = find(ct_labvec==n);
            testinds = randi(length(curcellinds),1,ceil(length(curcellinds)/k_));
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
        if savegroups_
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
    gmmstruct.std_err = std(error);
    gmmstruct.error = mean(error);
    gmmstruct.accuracy = mean(accuracy);
    gmmstruct.std_accuracy = std(accuracy);
    gmmstruct.gmmpost = mean(postprob);
    gmmstruct.std_gmmpost = std(postprob);
    gmmstruct.testgroups = group_mat;
    gmmstruct.truthgroups = test_mat;
else
    rng(0);
    centroids = zeros(size(C_ind_red,1),nts);
    C_testcells = [];
    test_labvec = [];
    C_traincells = [];
    for n = 1:nts
        curcellinds = find(ct_labvec==n);
        testinds = randi(length(curcellinds),1,ceil(length(curcellinds)/k_));
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
    if savegroups_
        group_mat = groups_ids;
        test_mat = test_labvec;
    else
        group_mat = [];
        test_mat = [];
    end
    accuracy = sum(test_labvec==groups_ids) / length(test_labvec);
    error = 1 - accuracy;
    gmmstruct.error = mean(error);
    gmmstruct.accuracy = mean(accuracy);
    gmmstruct.gmmpost = mean(post_probs);
    gmmstruct.testgroups = group_mat;
    gmmstruct.truthgroups = test_mat;
end

end