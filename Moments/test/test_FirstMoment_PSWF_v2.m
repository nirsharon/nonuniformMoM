% script name: "test_FirstMoment_PSWF_v2"
%

clear; clc;

% generate GT and gamma 
P = 3;
S = load('SO3_fifteen.mat');
SO_grid  = S.SO3;
gridSize = 15;
delta    = .999;

[A, B, vol, gamma] = createGTdata(P, gridSize, SO_grid, delta);
%K = size(gamma.coeff{1},2); %size(A{1}{1},1);


% preprocessing
L = max(gamma.band_idx_3d);  % the overall degree. Common to both 2D and 3D
P = size(B,1);               % expansion length of the distribution
M = max(gamma.ang_idx_2d)+1; % angular size of the moment

tic
[sign_mat, Gamma_mat] = preprocessing_mu1_coef(P, L);

[C_tensor] = preprocessing_mu2_coef(gamma, C_array, P);


% in future we should unit those two...
[C_array] = mu2_C_coefs_PSWF_v1_fixed(L, P, M);
[Gamma_mat, C_tensor] = preprocessing_mu2_coef(gamma, C_array, P);

tt = toc;
fprintf('DONE in about %d seconds \n', round(tt));

% first moment
fprintf('Mu1 calculation...');
tic
[m1] = FirstMoment_PSWF_v1_fixed(A, B, gamma);
tt1 = toc;

tic
[m1_new] = FirstMoment_PSWF_v2(A, B, Gamma_mat, sign_mat);
tt2 = toc;

norm(m1(:)-m1_new(:))

