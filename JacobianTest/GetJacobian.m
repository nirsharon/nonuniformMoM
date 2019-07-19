function [Jacobian_tensor] = GetJacobian(A, B, gamma, C_array)
% This function return the Jacobian of mu2 with respect to B and A
% A 5D-array J(m1, k1, m2, k2, :), last coordinate column corrsponds to
% vec_AB, the vectorial representation of A and B (together).

% parameters
K = size(gamma.coeff{1},2); %size(A{1}{1},1);
L = max(gamma.band_idx_3d);  % the overall degree. Common to both 2D and 3D
P = size(B,1);               % expansion length of the distribution
M = max(gamma.ang_idx_2d)+1; % angular size of the moment


% vectorize A and B
ind_Bvec = @(p,u,v) (2*p-1)*(2*p)*(2*p+1)/6 + u + (v-1)*(2*p+1);
vec_B  = FromCellArr2Vec({1},B);
vec_AB = A_B_to_VecAB(A, B, gamma.band_idx_3d);
size_vec_AB = length(vec_AB);
size_B = numel(vec_B)-1;
size_A = size_vec_AB - size_B;

% intialization
Jacobian_tensor = zeros(M, K, M, K, size_vec_AB);
A_coefs_tensor = zeros(M, K, M, K, size_A);
B_coefs_tensor = zeros(sizeM, sizeK, sizeM, sizeK, size_B);

% A part
for m1=0:(M-1)
    for m2=0:(M-1)
       current =  
       Jacobian_tensor(m1+1, k1_ind, m2+1, k2_ind, 1:size_A)  = current;
    end
end

% B part
B_coefs_tensor =  getting_B_coefs_PSWF(A, B, gamma, C_array);
Jacobian_tensor(:, :, :, :, (size_A+1):end) = B_coefs_tensor;

end

