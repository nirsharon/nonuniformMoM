% script name: generate_vol_dist 
%
% By BL, Apr 19

%% Create some probability (in alpha, beta, gamma). This probability will be replaced by a probability with finite expansion in Wigner-D
P  = 4;     % aribitrary, chosen to be small

% R1 = randn(3); [R1,~] = qr(R1,0);
% R2 = randn(3); [R2,~] = qr(R2,0);
R1 = eye(3);
R2 = eul2rotm([0,pi/2,0]);

if P==1
    B{1} = 1;
else
%     rho = @(a) exp(-5*norm(eul2rotm(a)-[0 1 0;0 0 1;1 0 0],'fro')^2)+ exp(-5*norm(eul2rotm(a)-[1 0 0;0 1 0;1 0 0],'fro')^2);
    rho = @(a) exp(-5*norm(eul2rotm(a)-R1,'fro')^2)+ 1*exp(-5*norm(eul2rotm(a)-R2,'fro')^2);
    B   = WignerD_transform(rho, P);  % the coefficients
    %"Project" to positive
    load('SO3_fifteen.mat');
    [AI,AE,Acon] = linear_cons_B2(numel(B)-1,SO3);
    [B,Breal]    = project_B_to_positive2(AI,AE,Acon,B);
    scl = 1/B{1};
    for i=1:numel(B)
        B{i} =  scl*B{i};
    end
    %     % test B
    %     rho_new = @(a) wignerD_expansion(B,a);
    %     B2 = WignerD_transform(rho_new, P);
    %     norm(blkdiag(B{:})-blkdiag(B2{:}))
end

%% Make distribution in-plane uniform
% for i=1:numel(B)
%     tmp = zeros(2*i-1);
%     tmp(:,i) = B{i}(:,i);
% %     tmp(:,i) = B{i}(i,:);
%     B{i} = tmp;
% end

%% Simulate kspace function
gridSize = 21;  % number of voxels in each dimenison. take odd.
% gridSize = 15;  % number of voxels in each dimenison. take odd.

% Sum of Gaussian
sigma = 200;
T     = (gridSize/5)*[0 0 0; 0.08 0.1 0; -0.1 0 0.1]';% ; 0.13 -0.2 0.1;0.1 0 -0.15]';
% T     = (gridSize/5)*[0 0 0; 0.12 0.12 0; 0 -0.12 0.1]';% ; 0.13 -0.2 0.1;0.1 0 -0.15]';
% T     = (gridSize/5)*[0 0 0; 0.3 0.3 0.3;-0.3 0 -0.3]';
% mu_x = 0.3;
% mu_y = 0.2;
% T     = (gridSize/2)*[0.3 0.2 0; 0.3 0.3 0.3; -0.3 0 -0.3]';
g     =  @(k) exp(1i*k'*T(:,1)).*exp(-pi^2/sigma*sum(k.*k)).'  + ...
    exp(1i*k'*T(:,2)).*exp(-pi^2/sigma*sum(k.*k)).'  + ...
    exp(1i*k'*T(:,3)).*exp(-pi^2/sigma*sum(k.*k)).' ;%+ ...
%          exp(1i*k'*T(:,4)).*exp(-pi^2/sigma*sum(k.*k)).'; % + ...
%          exp(1i*k'*T(:,5)).*exp(-pi^2/sigma*4*sum(k.*k)).';
r 