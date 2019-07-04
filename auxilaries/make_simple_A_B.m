% script name: "make_simple_A_B"
P  = 3;     % aribitrary, chosen to be small

rho = @(a) exp(-5*norm(eul2rotm(a)-[0 1 0;0 0 1;1 0 0],'fro')^2)+ exp(-5*norm(eul2rotm(a)-[1 0 0;0 0 1;1 0 0],'fro')^2);
B   = WignerD_transform(rho, P);  % the coefficients
%"Project" to positive
load('SO3_fifteen.mat');
[AI,AE,Acon] = linear_cons_B2(numel(B)-1,SO3);
[B,Breal]    = project_B_to_positive2(AI,AE,Acon,B);
scl = 1/B{1};
for i=1:numel(B)
    B{i} =  scl*B{i};
end


%% Simulate kspace function
gridSize = 21;  % number of voxels in each dimenison. take odd.

% Sum of Gaussian
sigma = 200;
T     = (gridSize/5)*[0 0 0; 0.08 0.1 0; -0.1 0 0.1]';% ; 0.13 -0.2 0.1;0.1 0 -0.15]';
g     =  @(k) exp(1i*k'*T(:,1)).*exp(-pi^2/sigma*sum(k.*k)).'  + ...
    exp(1i*k'*T(:,2)).*exp(-pi^2/sigma*sum(k.*k)).'  + ...
    exp(1i*k'*T(:,3)).*exp(-pi^2/sigma*sum(k.*k)).' ;%+ ...

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
