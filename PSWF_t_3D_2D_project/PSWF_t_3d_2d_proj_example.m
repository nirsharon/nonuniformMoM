clc;
clear variables
close all;

%% Define expansion parameters
beta = 1;       % Bandlimit ratio (between 0 and 1) - smaller values correspond for greater oversampling (choose beta=1 for Nyquist sampling, beta=0.5 for x2 oversampling, etc.)
delta = 0.99;   % Truncation parameter (between 0 and 1) - small values improve accuracy, large values use shorter expansions and improve noise spectrum. Approximation error should be small even for values near 1 if the volume is spatially localized. 
eps_p = 1e-3;   % Prescribed approximation accuracy for expanding 2D projection images using PSWFs.

%% Define test volume parameters
gridSize = 47;  % number of voxels in each dimenison - take odd
t = 0.04;       % Gaussian function standard deviation

%% Generate synthetic example: Non-centered Gaussian
L = floor(gridSize/2);
if mod(gridSize,2)==0
    x_1d = (-L:1:L-1)/L;   % - Even number of points
else
    x_1d = (-L:1:L)/L;   % - Odd number of points
end
[x_3d,y_3d,z_3d] = meshgrid(x_1d,x_1d,x_1d);
r = sqrt(x_3d.^2 + y_3d.^2 + z_3d.^2);
ball = (r <= 1);

mu_x = 0.3; 
mu_y = 0.2;
vol = exp(-((x_3d-mu_x).^2 + (y_3d-mu_y).^2 + z_3d.^2)/t);

%% Display central slices of volume function and Fourier magnitude
figure; imagesc(vol(:,:,L+1));
volFFT = fftshift(fftn(vol));
figure; imagesc(abs(volFFT(:,:,L+1)).^2);

%% Obtain 3D expansion coefficients
coeffs_3d = pswf_t_f_3d(vol, beta, delta);          

%% Reconstruct volume from coefficients
vol_hat = pswf_t_b_3d(coeffs_3d, gridSize, beta, delta);    

%% Compute 3D reconstruction error and display reconstructed slice
vol_err_calc_idx = 1;

vol_test = vol(:,:,:,vol_err_calc_idx);
vol_hat_test = vol_hat(:,:,:,vol_err_calc_idx);
figure; imagesc(vol_test(:,:,L+1));
figure; imagesc(vol_hat_test(:,:,L+1));

disp(['3D Gaussian reconstruction - relative squared error (inside the unit ball): ',num2str(mean(abs(vol_test(ball)-vol_hat_test(ball)).^2)/mean(abs(vol_test(ball)).^2))])

%% Approximate central slice in Fourier domain by DFT
% x_grid_3d = [x_3d(:),y_3d(:),z_3d(:)];
% [k_x_2d,k_y_2d] = meshgrid(x_1d);
% omega_grid = [k_x_2d(:) k_y_2d(:) zeros(numel(k_x_2d(:)),1)];
% DFT_mat = exp(-1i*pi*L*omega_grid*x_grid_3d.');
% proj_F_vec = DFT_mat*vol(:)/L^3;
% proj_F = reshape(proj_F_vec,gridSize,gridSize);

%% Approximate inverse Fourier transform by inverse DFT
% [x_2d,y_2d] = meshgrid(x_1d);
% x_grid_2d = [x_2d(:) y_2d(:)];
% iDFT_mat = exp(1i*pi*L*(x_grid_2d*x_grid_2d.'));
% proj_vec = iDFT_mat*proj_F(:)/4;
% proj = reshape(real(proj_vec),gridSize,gridSize);

%% Approximate central projection by simple numerical integration
% proj = sum(vol,3)/L;

%% Generate 2D projection analytically
[x_2d,y_2d] = meshgrid(x_1d,x_1d);
proj = sqrt(pi*t)*exp(-((x_2d-mu_x).^2 + (y_2d-mu_y).^2 )/t);

%% Compute mapping coefficients from 3D PSWF coefficients to 2D PSWF coefficients
c = beta*pi*L;  % Nyquist bandlimit
gamma = PSWF_2D_3D_T_mat(c, delta, eps_p);
bandMax = max(gamma.band_idx_3d);

%% Expand projection in 2D PSWFs
[coeffs_2d, PSWF_Nn_p_2d] = pswf_t_f_2d(proj, L, beta, eps_p, 1, []);

% Remark: Note that function pswf_t_f_2d outputs more coefficients than
% required accoring to the truncation rule by the gamma coefficients. This
% is because it uses a the more general truncation rule based on
% localization in the Disk.

%% Compare 2D PSWF coefficients with theory (i.e. according to projection from 3D)
a_mk_theory = [];
a_mk_empirical = [];
for m=0:max(gamma.ang_idx_2d)
    for k=0:nnz(gamma.ang_idx_2d==m)-1
        a_mk_theory_curr = 0;
        for l = abs(m):bandMax
            for s = 0:nnz(gamma.band_idx_3d==l)-1
                [coeff_lsmk] = get_Gamma_Coeff(l,s,m,k,gamma);
                a_mk_theory_curr = a_mk_theory_curr + coeff_lsmk * coeffs_3d{1}{l+1}(s+1,m+1);
            end
        end
        a_mk_theory = [a_mk_theory; a_mk_theory_curr];
    end
    a_mk_empirical_curr_m = coeffs_2d(PSWF_Nn_p_2d.ang_freq==m);
    a_mk_empirical = [a_mk_empirical; a_mk_empirical_curr_m(1:nnz(gamma.ang_idx_2d==m))];
end

e = norm(a_mk_empirical - a_mk_theory)^2/norm(a_mk_theory)^2;
disp(['Relative squared error between empirical 2D coefficients and theoretical 2D coefficients: ',num2str(e)]);
