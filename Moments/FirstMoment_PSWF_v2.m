function [ mu1 ] = FirstMoment_PSWF_v2(A, B, Gamma_mat, sign_mat, A_mat)
%
% Calculating the first moment in Prolates coordinates. 
%          ====    After fixing, April 19.   ===
% Input:
%    A       -- Cell array of the volume in 3D PSWF coordinates
%    B       -- Cell array {p,1}, each is a matrix of size (2p+1) X (2p+1)
%               (total number of coefficients = (2P-1)2P(2P+1)/6)
%    gamma   -- a struct that includes all info regarding the conversion between 
%                PSWF coefficients of 2D and 3D objects
% Output:
%     mu 1   -- a 2D array of the first moment, mu1 = mu1(m,k)
%
% NS, March 19
%------------------------------------------

% the maximal degree
max_deg = min(length(A{1})-1, size(B,1)-1);  

% total sizes of the moment
sizeM  = length(A{1});
sizeK  = size(Gamma_mat{1},1);

if nargin<5
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

% initialize mu1
mu1 = zeros(sizeM, sizeK);

% main loop
for m=0:max_deg     
    current_k = size(Gamma_mat{m+1},1);
    
    mu1(m+1,1:current_k) = Gamma_mat{m+1}*A_mat*sign_mat{m+1}*vec_B;
end

end

