function [mu2] = SecondMoment_PSWF_v2(A, B, gamma, C_tensor, Gamma_mat, A_mat)
%
%  This function evaluates the second moment in prolates coordinates, given 
%  the coefficients of the volume and distribution.
%
%          ====    After fixing, April 19.   ===
%
%  First "fast version", to use in simulations
%
% Input:
%    A       -- Cell array of the volume in 3D PSWF coordinates
%    B       -- Cell array {p,1}, each is a matrix of size (2p+1) X (2p+1)
%               (total number of coefficients = (2P-1)2P(2P+1)/6)
%    gamma   -- a struct that includes all info regarding the conversion between 
%                PSWF coefficients of 2D and 3D objects
%    C_array -- array of clebsch-gordon coefficients (preprocessing)
%
% Output:
%    mu2 -- the second moment in PSWF parameterization (m,k,m,k))
%
% NS, Feb 19
%------------------------------------------

% the overall degree. Common to both 2D and 3D
L = max(gamma.band_idx_3d);

% total sizes of the moment
sizeM  = max(gamma.ang_idx_2d)+1 ;
sizeK  = nnz(gamma.ang_idx_2d==0);
C_size = (L+1)^2; % old: (L^2+3*L+2)/2;

if nargin<6
    % move A to blkdiag
    A_mat = [];
    for m = 0:sizeM-1
        Al_s_n      = A{1}{m+1}(:,1:(m+1));
        current_blk = [ conj(Al_s_n(:,end:-1:2))*diag((-1).^(-m:1:-1).') , Al_s_n];
        A_mat = blkdiag(A_mat,current_blk);
    end
end

% vectorize B
vec_B = FromCellArr2Vec({1},B);

mu2_calc = zeros(sizeM, sizeM, sizeK, sizeK);
% main loop
for m1 = 0:sizeM-1
    for m2 = 0:sizeM-1
        
        
        %=========== THE MAIN PRODUCT IS DONE HERE ============
        Big_C = reshape(C_tensor{m1+1,m2+1}*vec_B, C_size, C_size);
        current_prod = Gamma_mat{m1+1}*A_mat*Big_C*(Gamma_mat{m2+1}*A_mat)';
        %=========== THE MAIN PRODUCT IS DONE HERE ============
        
        k1_ind = (1:nnz(gamma.ang_idx_2d==m1));
        k2_ind = (1:nnz(gamma.ang_idx_2d==m2));
        mu2_calc(m1+1,m2+1,k1_ind,k2_ind) = current_prod(k1_ind,k2_ind);
    end
end

% a final reshape
mu2 = zeros(sizeM, sizeK, sizeM, sizeK);
for  k1=1:sizeK
    mu2(:, k1, :, :) = mu2_calc(:, :, k1, :);
end

end


