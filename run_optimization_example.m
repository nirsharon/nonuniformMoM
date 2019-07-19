% script name: "run_optimization_example"
%
% This script to demonstrate the optimization code
%
% Be sure to have all supporting components in your path

saveit = 0; % to save the vols as mrc or not

% preliminaries
fprintf('\n  \n  \n  \n')
clear;

%% generating "ground truth"
fprintf('**** Starting Demonstration ****\n');
fprintf('\n  \n')
fprintf('Loading data..')
load('data_for_example')

% parameters
MAX_iter = 100;  % maximal iteration
gridSize = size(vol,1);   % volume resolution
beta     = c/(pi*radius);
delta    = .99;
[gamma_t,A_t] = gamma_truncate_2(gamma, A); % only for degrees with S(l)>=2l+1

% main parameters
L = max(gamma_t.band_idx_3d);  % the overall degree. Common to both 2D and 3D
P = size(B,1);               % expansion length of the distribution
M = max(gamma_t.ang_idx_2d)+1; % angular size of the moment
fprintf('Length of original expansion: %d and the truncated one is of: %d \n', length(A{1}), L)

% vectorized version
vec_AB_GT = A_B_to_VecAB(A_t, B, M-1);
vec_A_GT  = A_B_to_VecAB(A_t, [], M-1);
vec_B     = FromCellArr2Vec({1},B);

fprintf('Done  \n');

%% calculating moments
fprintf('Calculating moments: \n');

% preprocessing
tic
fprintf('Preprocessing..');
[sign_mat, Gamma_mat] = preprocessing_mu1_coef(P, L, gamma_t);
[C_tensor]            = preprocessing_mu2_coef_V2(gamma_t, P, M);
%[C_array] = mu2_C_coefs_PSWF_v1_fixed(L, P, M); % NO NEED. just for comparison with older version
tt = toc;
fprintf('DONE in about %d seconds \n', round(tt));

% first moment
[ m1_true ] = FirstMoment_PSWF_v2(A_t, B, Gamma_mat, sign_mat);
fprintf('First moment is DONE. \n');

% second moment
fprintf('Mu2 calculation...');
tic
[m2_true] = SecondMoment_PSWF_v2(A_t, B, gamma_t, C_tensor, Gamma_mat);
tt = toc;
fprintf('DONE in about %d seconds \n', round(tt));


%% setting the optimization
fprintf('Setting the optimization \n');

initial_guess = randn(size(vec_A_GT))*.001;
N = length(initial_guess);

m1 = m1_true;
m2 = m2_true;

%% running optimization
tic
manifold  = euclideancomplexfactory(N, 1);
problem.M = manifold;
weight    = numel(m1(:))/numel(m2(:));

problem.costgrad = @(x) LS_cost_w_grad_JUST_VOL_grad_fix(x, B, gamma_t, C_tensor, Gamma_mat , ...
    sign_mat, m1, m2, weight);

problem.hess = @(x, u)  getHessVal(x, B, gamma_t, C_tensor, Gamma_mat , ...
    sign_mat, m1, m2, weight, u);

options.maxinner          = 25;
options.tolgradnorm       = 1e-20;
options.maxiter           = MAX_iter;
[x, xcost, info, options] = trustregions(problem, initial_guess, options);


[A_est, ~] = VecAB_to_A_B(x, gamma_t);
toc();

t_opt = toc;


%% save the results

A_est_padded  = A;
A_padded      = A;
for j=1:length(A{1})
    if j<=length(A_est{1})
        A_est_padded{1}{j} = A_est{1}{j};
        A_padded{1}{j}     = A{1}{j};
    else
        A_est_padded{1}{j} = zeros(size(A_est_padded{1}{j}));
        A_padded{1}{j}     = zeros(size(A{1}{j}));
    end
end

%inverse prolates stransform
vol_hat  = pswf_t_b_3d(A_est_padded, gridSize, beta, delta);
vol_trnc = pswf_t_b_3d(A_padded,     gridSize, beta, delta);

%% if saveit~=0
if saveit
    nameit = 'simple_example';
    
    % print out volumes
    VisualVol(vol_hat,['est_',nameit])
    VisualVol(vol_trnc,['trunc_',nameit]);
    VisualVol(vol, ['origin_vol_',nameit]);
    [~, ~, volR3] = cryo_align_densities(vol, vol_hat, 3.35,1);
    VisualVol(volR3, ['rot_est_vol_',nameit]);
    % [resA,fighandle] = plotFSC(vol_hat, vol, .5, 3.35);
end
