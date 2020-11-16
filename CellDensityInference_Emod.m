function [B,res] = CellDensityInference_Emod(E_red,C_red,resflag)
% This function solves the nonnegative matrix inversion problem E = C*D,
% where E is the genes x voxels ISH data matrix, C is the genes x types
% scRNAseq data matrix, and B is the (to-be-determined) voxels x types
% inferred density matrix. The problem is solved one voxel at a time, as
% they can be treated independently. A note - in Mezias et al, 2020, the
% variable known as "B" is referred to as "D" and we retain that 
% nomenclature for backwards-compatibility purposes.

if nargin < 3
    resflag = '';
end

% Row (per-gene) normalization prior to matrix inversion. E should already
% be normalized but it is performed here for completeness.
N_g = size(C_red,1);
N_t = size(C_red,2); 
N_v = size(E_red,2);
for i = 1:N_g
    Crowmean = sum(C_red(i,:));
    C_red(i,:) = C_red(i,:)/Crowmean;
    Erowmean = sum(E_red(i,:));
    E_red(i,:) = E_red(i,:)/Erowmean;
end

% Matrix inversion using lsqnonneg
B = zeros(N_t,N_v);
res = zeros(N_t,N_v);
for i = 1:N_v
%     sprintf('Voxel %d/%d', i, N_v)
    if strcmp(resflag,'e')
        [curb,~,error] = lsqnonneg(C_red, E_red(:,i));
        B(:,i) = curb;
        res(:,i) = error;
    elseif strcmp(resflag,'e2')
        [curb,error] = lsqnonneg(C_red, E_red(:,i));
        B(:,i) = curb;
        res(:,i) = error;
    elseif strcmp(resflag,'r2')
        B(:,i) = lsqnonneg(C_red, E_red(:,i));
        xterm = C_red * B(:,i);
        yterm = E_red(:,i);
        res(:,i) = corr(xterm,yterm).^2;
    else
        options = optimoptions('lsqlin','Display','off');
        B(:,i) = lsqlin(C_red, E_red(:,i),[],[],[],[],zeros(N_t,1),Inf(N_t,1),[],options);
        res = 1;
    end
end
B = B.';
res = res.';
B = B .* res;
end