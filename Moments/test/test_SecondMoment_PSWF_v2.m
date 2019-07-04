% script name: "test_SecondMoment_PSWF_v2"
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
% in future we should unit those two...THE FUTURE IS HERE!
[C_array] = mu2_C_coefs_PSWF_v1_fixed(L, P, M); % NO NEED. just for comparison with older version
% [Gamma_mat, C_tensor] = preprocessing_mu2_coef(gamma, C_array, P); OLD!

[~, Gamma_mat] = preprocessing_mu1_coef(P, L, gamma);  % for gamma mat
[C_tensor]     = preprocessing_mu2_coef_V2(gamma, P, M);
tt = toc;
fprintf('DONE in about %d seconds \n', round(tt));

% second moment
fprintf('Mu2 calculation...');
tic
[m2] = SecondMoment_PSWF_v1_fixed(A, B, gamma, C_array);
tt1 = toc;

tic
[m2_prod] = SecondMoment_PSWF_v2(A, B, gamma, C_tensor, Gamma_mat);
tt2 = toc;

norm(m2(:)-m2_prod(:))

