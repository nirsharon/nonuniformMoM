% script name: "test_SecondMoment_PSWF_v1"

clear; clc;

% generate GT and gamma 
P        = 3;
S        = load('SO3_fifteen.mat');
SO_grid  = S.SO3;
gridSize = 15;
delta    = .99999999;

[A, B, vol, gamma] = createGTdata(P, gridSize, SO_grid, delta);

% preprocessing
L = max(gamma.band_idx_3d);  % the overall degree. Common to both 2D and 3D
P = size(B,1);               % expansion length of the distribution
M = max(gamma.ang_idx_2d)+1; % angular size of the moment

[C_array] = mu2_C_coefs_PSWF_v1_fixed(L, P, M);

tic
[mu2_v1] = SecondMoment_PSWF_v1_fixed(A, B, gamma, C_array);
%[mu2_v1] = SecondMoment_PSWF_v1(A, B, gamma, C_array);
toc();

% test
tic
[mu2] = SecondMoment_PSWF_naive_fixed(A, B, gamma);
toc();

norm(mu2(:)-mu2_v1(:))