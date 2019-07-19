function [Gamma_A] = GetGammaMat(gamma)
% This function calculates the gamma matrix that changes coordinates from
% 3D to 2D. Namely, 
%       Gamma_A(:,:,m) vec_A = vector of k associated with degree m

L = max(gamma.band_idx_3d);

% total sizes of the moment
sizeM = max(gamma.ang_idx_2d)+1 ;
sizeK = nnz(gamma.ang_idx_2d==0);
sizeA = sum(gamma.band_idx_3d)+ size(gamma.band_idx_3d,2);
% main loop - over the degrees (m is in 2d and L is in 3d)
Gamma_A = zeros(sizeK, sizeA ,sizeM);
for m1 = 0:sizeM-1
    start = 0;
        for l1 = m1:L
                gamma_s1_k1 = gamma.coeff{m1+1}(gamma.coeff_l_indices{m1+1}==l1,(1:nnz(gamma.ang_idx_2d==m1)));
               % Al1_s1_n1   = A{1}{l1+1}(:,1:(l1+1));
               % T_k1_n1     = gamma_s1_k1.'*Al1_s1_n1;
               if start == 0
                    start = sum(gamma.band_idx_3d(gamma.band_idx_3d<=(l1-1))+1)+1;
               else
                   start = ends +1;
               end
                ends = start + sum(gamma.band_idx_3d==l1) - 1;
                %ends  = sum(gamma.band_idx_3d(gamma.band_idx_3d<=l1)+1);
                s_ind_range = start:ends;
                k_ind = (1:nnz(gamma.ang_idx_2d==m1));
                %start = sum(gamma.band_idx_3d(gamma.band_idx_3d<=(l1-1))+1)+1;

               % Gamma_A(k_ind,s_ind_range, m1+1) = Gamma_A(k_ind,s_ind_range, m1+1) + repmat(gamma_s1_k1.',1,l1+1);
                Gamma_A(k_ind,s_ind_range, m1+1) = gamma_s1_k1.';
        end
end

end

