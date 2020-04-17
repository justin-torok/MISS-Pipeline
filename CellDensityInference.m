function [B] = CellDensityInference(E_red,C_red)
% This function solves the nonnegative matrix inversion problem E = C*D,
% where E is the genes x voxels ISH data matrix, C is the genes x types
% scRNAseq data matrix, and B is the (to-be-determined) voxels x types
% inferred density matrix. The problem is solved one voxel at a time, as
% they can be treated independently. A note - in Mezias et al, 2020, the
% variable known as "B" is referred to as "D" and we retain that 
% nomenclature for backwards-compatibility purposes.

% Row (per-gene) normalization prior to matrix inversion. E should already
% be normalized but it is performed here for completeness.
N_g = size(C_red,1);
N_t = size(C_red,2); 
N_v = size(E_red,2);
for i = 1:N_g
    Crowsum = sum(C_red(i,:));
    C_red(i,:) = C_red(i,:)/Crowsum;
%     Erowsum = sum(E_red(i,:));
%     E_red(i,:) = E_red(i,:)/Erowsum;
end

% Matrix inversion using lsqnonneg
B = zeros(N_t,N_v);
for i = 1:N_v
%     sprintf('Voxel %d/%d', i, N_v)
    B(:,i) = lsqnonneg(C_red, E_red(:,i));
end
B = B.';
end