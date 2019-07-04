% script name: "test_mu2_convergence_v6"
%
% We test the moments code using simulated volume in k-space
%
% Note: the most expensive part is calculating the second moment, with a new code
% we should repeat the test for larger P, for example.
%
% NS, 2019

clear all; close all;

%% Create some probability (in alpha, beta, gamma). This probability will be replaced by a probability with finite expansion in Wigner-D
P  = 2;     % aribitrary, chosen to be small

B{1} = 1;
val1 = rand*.01*1i;
B2(1,1) = val1; B2(3,1) = val1;
B2(1,3) = - B2(1,1); B2(3,3) = -B2(3,1);
val2 = rand*.01;
B2(1,2) = val2; B2(3,2) = -val2; 
B2(2,1) = -1i*val2; B2(2,3) = B2(2,1); 
B{2} = B2;
%     % "Project" to positive
%     load('SO3_fifteen.mat');
%     [AI,AE,Acon] = linear_cons_B2(numel(B)-1,SO3);
%     [BB,Breal]    = project_B_to_positive2(AI,AE,Acon,B);
%     scl = 1/B{1};
%     for i=1:numel(B)
%         B{i} =  scl*B{i};
%     end
%     %     % test B
%     %     rho_new = @(a) wignerD_expansion(B,a);
%     %     B2 = WignerD_transform(rho_new, P);
%     %     norm(blkdiag(B{:})-blkdiag(B2{:}))

%% Simulate kspace function
gridSize = 21;  % number of voxels in each dimenison. take odd.

% Sum of Gaussian
sigma = 200;
T     = (gridSize/5)*[0 0 0; 0.08 0.1 0; -0.1 0 0.1]';% ; 0.13 -0.2 0.1;0.1 0 -0.15]';
g     =  @(k) exp(1i*k'*T(:,1)).*exp(-pi^2/sigma*sum(k.*k)).'  + ...
    exp(1i*k'*T(:,2)).*exp(-pi^2/sigma*sum(k.*k)).'  + ...
    exp(1i*k'*T(:,3)).*exp(-pi^2/sigma*sum(k.*k)).' ;%+ ...
%          exp(1i*k'*T(:,4)).*exp(-pi^2/sigma*sum(k.*k)).'; % + ...
%          exp(1i*k'*T(:,5)).*exp(-pi^2/sigma*4*sum(k.*k)).';

%% Vol on kspace over cartesian grid

radius   = floor(gridSize/2);
if mod(gridSize,2)==0
    x_1d = (-radius:1:radius-1);%/radius;   % - Even number of points
else
    x_1d = (-radius:1:radius);%/radius;   % - Odd number of points
end

% the 3d grid
[x_3d,y_3d,z_3d] = meshgrid(x_1d,x_1d,x_1d);
vec_3d_grid      = [x_3d(:),y_3d(:),z_3d(:)].'*pi/2;

% parameterization of the 2d slices
x_2d = x_3d(:,:,1);
y_2d = y_3d(:,:,1);

% evaluate the function
volf = g(vec_3d_grid);

% back to real domain
volk = reshape(volf, gridSize, gridSize, gridSize);
vol  = real(fftshift(ifftn(ifftshift(volk)))); % in real domain
vol  = vol*floor(size(x_2d,1)/2);   % Nir Sunday ?

%% calculate 3D expansion and gamma coefficients
beta  = 1;       % Bandlimit ratio
delta = 0.99;    % Truncation parameter
eps_p = 1e-3;    % Prescribed accuracy

fprintf('Calculating 3D coefficients...');
A = pswf_t_f_3d(vol, beta, delta);
fprintf('DONE \n');

fprintf('Calculating gamma coefficients...');
c     = beta*pi*radius;              % Nyquist bandlimit
gamma = PSWF_2D_3D_T_mat(c, delta, eps_p);
fprintf('DONE \n');

%% checking reconstruct volume from coefficients
vol_hat = pswf_t_b_3d(A, gridSize, beta, delta);
% norm(vol(:)-vol_hat(:)) % norm(vol(:)-vol_hat(:))/norm(vol(:))

r    = sqrt(x_3d.^2 + y_3d.^2 + z_3d.^2);
ball = (r <= max(x_3d(:)));
% vol_supp = zeros(size(vol)); vol_supp(ball) = vol(ball); % figure; vol3d('cdata',real(vol_supp))

figure; vol3d('cdata',vol);     title('Original Volume')
figure; vol3d('cdata',vol_hat); title('Reconstructed Volume')
v1 = vol(ball);
v2 = vol_hat(ball);
disp(['3D vol reconstruction - relative squared error (inside the ball): ',num2str( (norm(v1(:)-v2(:)).^2)/norm(v1(:)).^2 ) ])

% % checking projection
% k_cord = [x_2d(:), y_2d(:), zeros(size(y_2d(:))) ].'*pi/2;
% read the slice in k-space
% current_im = g(k_cord);
% current_im = reshape(current_im, size(x_2d));
% proj       = real(fftshift(ifft2(ifftshift(current_im))));
% figure; imagesc(sum(vol,3)); figure; imagesc(proj);
% norm(sum(vol,3)-proj)
% norm(sum(vol,3)/floor(size(x_2d,1)/2)-proj)
% AA = randn(3);
% if det(AA)<0
%     ind = randi(3);
%     AA(:,ind) = -AA(:,ind);
% end
% [U,~,V]    = svd(AA); RR = U*V';
% R_im    = g(RR*k_cord);
% R_im    = reshape(R_im, size(x_2d));
% R_proj  = real(fftshift(ifft2(ifftshift(R_im))));
% rotating the vol
% axang = vrrotmat2vec(RR);
% R_vol = imrotate3(vol,axang(4),axang(1:3),'cubic','crop');
%
% figure; vol3d('cdata',R_vol);     title('Rotated Volume')
% figure; vol3d('cdata',vol);     title('Original Volume')
%
% figure; imagesc(sum(R_vol,3)); figure; imagesc(R_proj);
% norm(sum(R_vol,3)-R_proj)


%% calculating moments
fprintf('Calculating moments: \n');

% first moment
[ m1_true ] = FirstMoment_PSWF_v1_fixed(A, B, gamma);  % old code: [ m1_true ] = FirstMoment_PSWF_naive(A, B, gamma);
fprintf('First moment is DONE. \n');

% preprocessing
L = max(gamma.band_idx_3d);  % the overall degree. Common to both 2D and 3D
P = size(B,1);               % expansion length of the distribution
M = max(gamma.ang_idx_2d)+1; % angular size of the moment
fprintf('Run preprocessing...');
tic
[C_array] = mu2_C_coefs_PSWF_v1_fixed(L, P, M);
tt = toc;
fprintf('DONE in about %d seconds \n', round(tt));

% second moment
fprintf('Mu2 calculation...');
tic
[m2_true] = SecondMoment_PSWF_v1_fixed(A, B, gamma, C_array);
tt = toc;
fprintf('DONE in about %d seconds \n', round(tt));


%% prepare the comparison data
total_N      = 10000;
repeat_trial = 15;

fprintf('Preparing projections: \n');

% sampling SO(3))
fprintf('Uniform sampling...')
R = generate_SO3_uniform_array_v2(total_N);
fprintf('DONE \n');

% projection simulating
fprintf('projection simulating...');
[projs, weight] = make_projs_v2(g, R, B, x_2d, y_2d);
weight = real(weight);
weight(weight<0) = 0;
fprintf('DONE \n');

% overall energy
norm_m1_square = norm(m1_true(:))^2;
norm_m2_square = norm(m2_true(:))^2;

% comparison points
number_of_check_points = 50;
start_N                = 5;
check_points           = floor(linspace(start_N, total_N, number_of_check_points));

%% moving to prolates coordinates -- define expansion parameters
fprintf('Moving generated projections to PSWF coordinates...');
beta  = 1;      % Bandlimit ratio
eps_p = 1e-3;   % Prescribed accuracy
[im_coef, p_struct] = pswf_t_f_2d(projs, floor(size(projs,1)/2), beta, eps_p, 1, []);

% accuracy check for the 2d
%images = pswf_t_b_2d( im_coef, p_struct, 1 );
%norm(images(:,:,1)-projs(:,:,1))/norm(projs(:,:,1))
%norm(images(:)-projs(:))/norm(projs(:))

proj_PSWF = zeros(max(gamma.ang_idx_2d)+1, nnz(gamma.ang_idx_2d==0), size(im_coef,2));
for m=0:max(gamma.ang_idx_2d)
    for k=0:nnz(gamma.ang_idx_2d==m)-1
        coeff2Dcurr_m = im_coef(p_struct.ang_freq==m,:);
        proj_PSWF(m+1,k+1,:) = coeff2Dcurr_m(k+1,:);
    end
end
fprintf('...DONE \n')

%% average to attain the empirical moments
fprintf('Averaging...');
% initiliaze
relative_m1 = zeros(size(check_points));
relative_m2 = zeros(size(check_points));
absolute_m1 = zeros(size(check_points));
absolute_m2 = zeros(size(check_points));
m1_hat = 0;
m2_hat = 0;
for rt=1:repeat_trial
    % initialize inner variables
    check_counter = 1;
    inner_relative_m1 = zeros(size(check_points));
    inner_relative_m2 = zeros(size(check_points));
    inner_absolute_m1 = zeros(size(check_points));
    inner_absolute_m2 = zeros(size(check_points));
    perm_list = randperm(size(R,3));
    m1_hat = 0;
    m2_hat = 0;
    for i=1:total_N
        % averaging
        m1_hat = m1_hat + proj_PSWF(:,:,perm_list(i))*weight(perm_list(i));
        current_vec_proj = proj_PSWF(:,:,perm_list(i));
        m2_hat = m2_hat + current_vec_proj(:)*current_vec_proj(:)'*weight(perm_list(i));
        
        % check point control
        if ismember(i,check_points)
            
            % first moment
            current_m1_check     = m1_hat/sum(weight(perm_list(1:i)));
            inner_absolute_m1(check_counter) = norm(m1_true(:)-current_m1_check(:))^2;
            inner_relative_m1(check_counter) = inner_absolute_m1(check_counter)/norm_m1_square;
            
            % second moment initialization
            current_m2_check = zeros(max(gamma.ang_idx_2d)+1, nnz(gamma.ang_idx_2d==0),max(gamma.ang_idx_2d)+1, nnz(gamma.ang_idx_2d==0));
            
            %reshape the second moment and summarize
            for m=0:max(gamma.ang_idx_2d)
                for k=0:nnz(gamma.ang_idx_2d==0)-1
                    current_m2_check(m+1,k+1,:,:) = reshape(m2_hat(m+1 + (k)*(max(gamma.ang_idx_2d)+1),:),...
                        max(gamma.ang_idx_2d)+1, nnz(gamma.ang_idx_2d==0) );
                end
            end
            current_m2_check = current_m2_check/sum(weight(perm_list(1:i)));
            inner_absolute_m2(check_counter) = norm(m2_true(:)-current_m2_check(:))^2;
            inner_relative_m2(check_counter) = inner_absolute_m2(check_counter)/norm_m2_square;
            check_counter = check_counter + 1;        
        end 
    end
    % adding current trial results
    relative_m1 = relative_m1 + (1/repeat_trial)*inner_relative_m1;
    absolute_m1 = absolute_m1 + (1/repeat_trial)*inner_absolute_m1;
    relative_m2 = relative_m2 + (1/repeat_trial)*inner_relative_m2;
    absolute_m2 = absolute_m2 + (1/repeat_trial)*inner_absolute_m2;    
    fprintf('*');
end

fprintf('...DONE \n');

%% visualize the result
figure; plot(check_points/1000,relative_m1,'LineWidth',2.5); legend('\mu_1')
xlabel('Number of samples (in thousands)'); ylabel('Square relative error'); set(gca,'FontSize',22)

figure; plot(check_points/1000,relative_m2,'r','LineWidth',3); legend('\mu_2');
xlabel('Number of samples (in thousands)'); ylabel('Square relative error'); set(gca,'FontSize',22)

figure; plot(check_points/1000,absolute_m1,'LineWidth',2.5); legend('\mu_1')
xlabel('Number of samples (in thousands)'); ylabel('Square error'); set(gca,'FontSize',22)

figure; plot(check_points/1000,absolute_m2,'r','LineWidth',3); legend('\mu_2');
xlabel('Number of samples (in thousands)'); ylabel('Square error'); set(gca,'FontSize',22)
% tt = toc;
% fprintf('Total running time is about %d seconds \n', round(tt));

%% loglog
figure; loglog(check_points/1000,relative_m1,'LineWidth',2.5); legend('\mu_1')
xlabel('Number of samples (in thousands)'); ylabel('Square relative error'); set(gca,'FontSize',22)

figure; loglog(check_points/1000,relative_m2,'r','LineWidth',3); legend('\mu_2');
xlabel('Number of samples (in thousands)'); ylabel('Square relative error'); set(gca,'FontSize',22)
