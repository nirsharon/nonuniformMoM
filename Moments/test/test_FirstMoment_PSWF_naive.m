% script name: "test_FirstMoment_PSWF_naive"

clear; clc;

% generate GT and gamma 
P = 3;
S = load('SO3_fifteen.mat');
SO_grid  = S.SO3;
gridSize = 15;
delta    = .9999;

[A, B, vol, gamma] = createGTdata(P, gridSize, SO_grid, delta);
%K = size(gamma.coeff{1},2); %size(A{1}{1},1);

% test
tic
[mu1] = FirstMoment_PSWF_naive_fixed2(A, B, gamma);
toc();

tic
[mu1_v1] = FirstMoment_PSWF_v1_fixed(A, B, gamma);
toc();

norm(mu1(:)-mu1_v1(:))