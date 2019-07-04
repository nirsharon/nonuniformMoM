function [A, B, vol, gamma] = createGTdata(P, gridSize, SO_grid, delta, beta)
% Generating ground truth vol and dist for simulations
%
% Input:
%   gridSize - odd number of voxels in each dimenison of the volume

assert(mod(gridSize,2)==1, 'Wrong grid size, set to odd number');

if nargin<3
    S = load('SO3_fifteen.mat');
    SO_grid = S.SO3;
end

if nargin<4
    delta = 0.999;   % Truncation parameter
    beta  = 1;       % Bandlimit ratio
end

%% volume side
if gridSize>1
radius = floor(gridSize/2);
if mod(gridSize,2)==0
    x_1d = (-radius:1:radius-1)/radius;   % - Even number of points
else
    x_1d = (-radius:1:radius)/radius;   % - Odd number of points
end
[x_3d,y_3d,z_3d] = meshgrid(x_1d,x_1d,x_1d);
r = sqrt(x_3d.^2 + y_3d.^2 + z_3d.^2);
ball = (r <= 1);

mu_x = 0.3; 
mu_y = 0.2;
mu_z = 0.3; 

t   = 0.04;       % Gaussian function standard deviation
vol = exp(-((x_3d-mu_x).^2 + (y_3d-mu_y).^2 + (z_3d-mu_z).^2)/t) + ...
       exp(-((x_3d+mu_x-.05).^2 + (y_3d+mu_y+.1).^2 + z_3d.^2)/t)+ ...
       exp(-((x_3d+mu_x-.2).^2 + (y_3d+mu_y-.1).^2 + (z_3d.^2-.3).^2)/t)+ ...
       exp(-((x_3d+mu_x+.05).^2 + (y_3d+mu_y-.6).^2 + (z_3d+0.4).^2)/t);

% Obtain 3D expansion coefficients
eps_p = 5e-4;    % Prescribed accuracy . 

A     = pswf_t_f_3d(vol, beta, delta);
c     = beta*pi*radius;              % Nyquist bandlimit
gamma = PSWF_2D_3D_T_mat(c, delta, eps_p);

% Reconstruct volume from coefficients, if needed: vol_hat = pswf_t_b_3d(A, gridSize, beta, delta);    
else
    vol = [];
    gamma = [];
    A = [];
end

%% distribution side
B = cell(P,1);
pert_mag = .1;
B{1} = 1;
for j = 2:P
    B{j} = pert_mag*(rand(2*j-1)*.1 + rand(2*j-1 )*.1*1i)/(2*j-1) ;
end
[AI, AE, Acon] = linear_cons_B2(P-1,SO_grid);
[B, ~]         = project_B_to_positive2(AI, AE, Acon, B);
scl   = 1/B{1};
for i = 1:numel(B)
    B{i} =  scl*B{i};
end

end

