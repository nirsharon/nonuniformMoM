function [sign_mat, Gamma_mat] = preprocessing_mu1_coef(P, L, gamma)
%
% This function generates the sign matrix (array) which is the coefficient of the
% distribution in the first moment. After fixing, April 19.
%
% changing coordinates and sign matrix together!
%
% Input:
% P - length of distribution expansion, size(B,1)
% L - max degree of volume expansion, max(gamma.band_idx_3d);
%
% Output:
% sign_mat - cell array of matrices, according to the degree
% Gamma_mat - concat of gamma to form the changing corrdinates matrix
%
% NS, March 19

% the max degree
max_deg = min(L, P-1);

% size of the matrix
C_size = (L+1)^2; % old: (L^2+3*L+2)/2;

% length of the vectorized distribution, length(FromCellArr2Vec({1},B))
vec_B_len  = (2*P-1)*(2*P)*(2*P+1)/6;

% index functions
ind_Bvec  = @(p,u,v) (2*p-1)*(2*p)*(2*p+1)/6 + u + (v-1)*(2*p+1);
l_max_ind = @(l) (l+1)^2; %

% % initialization
% sign_mat2 = cell(max_deg+1,1);
% 
% % main loop
% for m = 0:max_deg
%     sign_mat2{m+1} = zeros(C_size, vec_B_len);
%     for l = m:max_deg
%         l_ind_rang = (l_max_ind(l-1)+1):(l_max_ind(l));
%         b_ind  = ind_Bvec(l, l+1+ (l:(-1):(-l)), l+1 -m);
%         sm_ind = sub2ind([C_size, vec_B_len], l_ind_rang, b_ind);
%         sign_mat2{m+1}(sm_ind) = (-1)^m*(-1).^((-l:l))/(2*l+1);
%     end
% end

% initialization
sign_mat  = cell(max_deg+1,1);
Gamma_mat = cell(L+1,1);

% main loop
for m = 0:L
    if m<=max_deg
        sign_mat{m+1} = zeros(C_size, vec_B_len);
    end
    for l = m:L
        if and((m<=max_deg),(l<=max_deg))
            l_ind_rang = (l_max_ind(l-1)+1):(l_max_ind(l));
            b_ind  = ind_Bvec(l, l+1+ (l:(-1):(-l)), l+1 -m);
            sm_ind = sub2ind([C_size, vec_B_len], l_ind_rang, b_ind);
            sign_mat{m+1}(sm_ind) = (-1)^m*(-1).^((-l:l))/(2*l+1);
        end
        gamma_s_k  = gamma.coeff{m+1}(gamma.coeff_l_indices{m+1}==l,(1:nnz(gamma.ang_idx_2d==m)));
        zero_count = find(gamma.coeff_l_indices{1}==m,1)-1;
        s_ind      = [zeros(1,zero_count)==1,gamma.coeff_l_indices{m+1}==l];
        Gamma_mat{m+1}(1:size(gamma_s_k,2),s_ind) = gamma_s_k.';
    end
end

end

