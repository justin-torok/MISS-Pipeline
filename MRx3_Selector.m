function geneinds = MRx3_Selector(C_raw,E,n,lambda)
% This function calculates the best n genes (rows) of matrix C according to
% the MRx3 criterion

% Create a column-normalized (i.e. by cell type) C matrix for the purpose
% of calculating the approximate F-statistic and pairwise correlations for
% the MRx3 criterion
ctmean = mean(C_raw,1);
ctmean = repmat(ctmean,size(C_raw,1),1);
ctnorm = C_raw ./ ctmean;
ctnorm(ctnorm<(0.1*std(nonzeros(ctnorm)))) = 0; % filter out excessively small signal
C_cellnorm = ctnorm;

% Create a row-normalized (i.e. by gene) C matrix for the purpose of
% calculating the rank-1 approximation to the projection error for the MRx3
% criterion
gsum = sum(C_raw,2);
gsum = repmat(gsum,1,size(C_raw,2));
gnorm = C_raw ./ gsum;
gnorm(isnan(gnorm)) = 0;
C_genenorm = gnorm;

Fcalc = @(g) sum((g - mean(g)).^2)/(length(g) - 1);
Rmat = corr(C_cellnorm.');
k_g = 10; % number of genes to add at a time (saves computational time)
for i = 1:size(E,2)
    curgene = E(:,i);
    cursum = sum(curgene);
    E(:,i) = curgene/cursum;
end
E = E.';

% mRMR Initialization
geneinds = [];
remaininds = 1:size(C_cellnorm,1);
Vinit = zeros(1,size(C_cellnorm,1));
for i = 1:size(C_cellnorm,1)
    g_ct = C_cellnorm(i,:);
    Vinit(i) = Fcalc(g_ct);
end
[Vmax,maxind] = max(Vinit);
geneinds = [geneinds,remaininds(maxind)];
remaininds(maxind) = [];
S_ct = C_cellnorm(geneinds,:);
S_g = C_genenorm(geneinds,:);
T_ct = C_cellnorm(remaininds,:);
T_g = C_genenorm(remaininds,:);

% mRMR Propagation/Termination
while length(geneinds) < n
%     fprintf('gene %d/%d\n', length(geneinds)+1,n)
    szT = size(T_ct,1);
    Vtest = zeros(1,szT);
    Fcalcvec = Vtest; wcalc = Vtest; errvec = Vtest;
    if length(geneinds) >= 30
        pinv_S = pinv(S_g);
        E_g = E(geneinds,:);
        D_est = pinv_S * E_g; 
    end
    for i = 1:szT
        g_ct = T_ct(i,:);
        g_g = T_g(i,:);
        if all(g_g == 0)
            Vtest(i) = -1000;
        else
            ind = remaininds(i);
            Rs = abs(Rmat(ind,geneinds));
            wc = sum(Rs)/length(geneinds);
            if length(geneinds) < 30 % Run normal mRMR for the first 30 genes
                Vtest(i) = Fcalc(g_ct) / wc;
            else
                E_i = E(ind,:);
                err1_i = sum((E_i - g_g*D_est).^2);
                err2_i = sum((g_g*pinv_S).^2)*dot(E_i,E_i);
                Fcalcvec(i) = Fcalc(g_ct); wcalc(i) = wc; errvec(i) = err1_i + err2_i;
                Vtest(i) = Fcalc(g_ct) / (wc + lambda*(err1_i + err2_i));
            end
        end
    end

    if length(geneinds) < 30
        [Vm,maxind] = max(Vtest);
        Vmax = [Vmax,Vm];
        geneinds = [geneinds,remaininds(maxind)];
        remaininds(maxind) = [];
    else
        Vtest(isnan(Vtest)) = 0;
        [Vm,sortind] = sort(Vtest,'descend');
        if (n - length(geneinds)) < k_g
            k_g = n - length(geneinds);
        end
        sortind = sortind(1:k_g);
        Vmax = [Vmax,Vm(1:k_g)];
        geneinds = [geneinds,remaininds(sortind)];
        remaininds(sortind) = [];
    end
    S_ct = C_cellnorm(geneinds,:);
    S_g = C_genenorm(geneinds,:);
    T_ct = C_cellnorm(remaininds,:);
    T_g = C_genenorm(remaininds,:);
end

end