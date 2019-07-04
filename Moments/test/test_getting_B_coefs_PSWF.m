% script name: "test_getting_B_coefs_PSWF"
%

clear; clc;

% generate GT and gamma 
P        = 4;
S        = load('SO3_fifteen.mat');
SO_grid  = S.SO3;
gridSize = 15;
delta    = .99999;

[A, B, vol, gamma] = createGTdata(P, gridSize, SO_grid, delta);

% preprocessing
L = max(gamma.band_idx_3d);  % the overall degree. Common to both 2D and 3D
P = size(B,1);               % expansion length of the distribution
M = max(gamma.ang_idx_2d)+1; % angular size of the moment

tic
fprintf('Preprocessing..');
[sign_mat, Gamma_mat] = preprocessing_mu1_coef(P, L, gamma);
[C_tensor]            = preprocessing_mu2_coef_V2(gamma, P, M);
[C_array] = mu2_C_coefs_PSWF_v1_fixed(L, P, M); % NO NEED. just for comparison with older version

tt = toc;
fprintf('DONE in about %d seconds \n', round(tt));

% second moment
tic
[m2_prod] = SecondMoment_PSWF_v2(A, B, gamma, C_tensor, Gamma_mat);
tt2 = toc;

% now with the B_coef
B_coefs_tensor = getting_B_coefs_PSWF_V2(A, B, gamma, C_array);

% summing up
vec_B = FromCellArr2Vec({1},B);
m2 = zeros(size(m2_prod));
for j=1:length(vec_B)
    m2 = m2 + B_coefs_tensor(:,:,:,:,j)*vec_B(j);
end

norm(m2(:)-m2_prod(:))

