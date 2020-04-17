function geneinds = mRMR_Selector(C,n,method)
% This function calculates the best n genes (rows) of matrix C according to
% an mRMR criterion

Fcalc = @(g) sum((g - mean(g)).^2)/(length(g) - 1);
Rmat = corr(C.');

% mRMR Initialization
geneinds = [];
remaininds = 1:size(C,1);
Vinit = zeros(1,size(C,1));
for i = 1:size(C,1)
    g = C(i,:);
    Vinit(i) = Fcalc(g);
end
[Vmax,maxind] = max(Vinit);
geneinds = [geneinds,remaininds(maxind)];
remaininds(maxind) = [];
% S = C(geneinds,:);
T = C(remaininds,:);

% mRMR Propagation/Termination
while length(geneinds) < n
%     sprintf('gene %d/%d', length(geneinds)+1,n)
    szT = size(T,1);
    Vtest = zeros(1,szT);
    for i = 1:szT
        g = T(i,:);
        ind = remaininds(i);
        Rs = abs(Rmat(ind,geneinds));
        wc = sum(Rs)/length(geneinds);
        if strcmp(method,'Diff')
            Vtest(i) = Fcalc(g) - wc;
        elseif strcmp(method,'Quo')
            Vtest(i) = Fcalc(g) / wc;
        end
    end
    [Vm,maxind] = max(Vtest);
    Vmax = [Vmax,Vm];
    geneinds = [geneinds,remaininds(maxind)];
    remaininds(maxind) = [];
%     S = C(geneinds,:);
    T = C(remaininds,:);
end

end