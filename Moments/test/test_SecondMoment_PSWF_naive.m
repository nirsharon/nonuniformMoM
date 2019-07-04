% script name: "test_SecondMoment_PSWF_naive"

clear; clc;

% load coeffs and gamma 
if exist('small_example.mat')
    load('small_example.mat');
else
    prepare_coef_and_gamma_example;
    clear;
    load('small_example.mat');
end

% random distribution
P = 2;
B    = cell(P,1);
B{1} = 1;
for j=2:P
    B{j} = randn(2*j-1,2*j-1);
end

K = size(gamma.coeff{1},2); %size(A{1}{1},1);

% test
tic
[mu2] = SecondMoment_PSWF_naive_fixed(A, B, gamma);
toc();

