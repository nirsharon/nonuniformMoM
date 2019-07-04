function [vol_hat] = pswf_t_b_3d(coeff, gridSize, beta, delta)
%   This is the backward 3D PSWF transform - from PSWF expansion coefficients to volumes sampled on the Cartesian grid.
%   The 3D PSWFs are given by \psi_{N,n,m}(r,\theta,\phi) = R_{N,n}(r) * Y_N^m(\theta,\phi),
%   where R_{N,n}(r) are radial functions and Y_N^m(\theta,\phi) are the spherical harmonics.
%   Input:  coeff:     PSWF expansion coefficients, provided as a cell
%                       array enumrated by the index of the volume function. Within each cell value, there is another cell array enumerated by the angular index N (equall to \ell in the classical spherical harmonics notation), where each cell contains a
%                       matrix of expansion coefficients such that the first dimension (rows)
%                       correspond to the radial index n, and the second dimension (columns)
%                       correspond to the spherical harmonic index m (which is between 0 and N since the volume is assumed to be real-valued).
%           gridSize:   Number of voxels for each dimension of the volume.
%           beta:       Bandlimit ratio relative to the Nyquist rate, between 0 and 1.
%           delta:      Truncation parameter, between 0 and 1, which controls the length of the
%                       expansion and the approximation error. Smaller values (closer to zero) guarantee smaller errors, yet longer expansions, and vice-versa. 
%                       Note: Due to numerical considerations, do not exceed [1E-6,1-1E-6]
%   Output: vol:        Reconstructed volume functions, provided as a 4D array. The fourth dimension enumerates over the different volumes, where the
%                       first three dimensions correspond to an equally-spacsed grid of gridSize x gridSize x gridSize points.
%
%
%   Boris Landa, 04.12.2018
%
%% Compute expansion coeffs for every angular frequency
L = floor(gridSize/2);
if mod(gridSize,2)==0
    x_1d = -L:1:L-1;   % - Even number of points
else
    x_1d = -L:1:L;   % - Odd number of points
end
[x_3d,y_3d,z_3d] = meshgrid(x_1d,x_1d,x_1d);
r = sqrt(x_3d(:).^2 + y_3d(:).^2 + z_3d(:).^2);
ball = (r <= L);
x_3d_b = x_3d(ball); y_3d_b = y_3d(ball); z_3d_b = z_3d(ball);

r_b = r(ball);
[r_b_unique,~,r_map] = unique(r_b);
theta_b = (acos(z_3d_b(:)./r_b));
phi_b =  (angle(x_3d_b(:)/L + 1i*y_3d_b(:)/L));

nVols = numel(coeff);
vol_mat_hat = zeros(nnz(ball),nVols);

c = beta * pi*L;

N = 0;
n = 2*L;
n_order_length_vec = [];
while(1)
    [PSWF_R,alpha] = PSWF_3D(N,n,c,r_b_unique/L,eps);
    PSWF_R = PSWF_R(r_map,:);
    
    lambda = (c/(2*pi))^3 * abs(alpha).^2;
    n_end = find(lambda<delta,1,'first');
    
    if (n_end < 2)
        break;
    end
        
    if (~isempty(n_end))
%         PSWF_R = PSWF_R*diag(abs(alpha).^2)*(beta/2)^3;
        PSWF_R = PSWF_R(:,1:n_end-1);
        Phase_part = mySph_v2(N,0:N,theta_b,phi_b);     % Complex valued basis functions
        for i=1:nVols
            vol_mat_hat(:,i) = vol_mat_hat(:,i) + real((PSWF_R*coeff{i}{N+1}(:,1)).*Phase_part(:,1));
            vol_mat_hat(:,i) = vol_mat_hat(:,i) + 2*sum(real((PSWF_R*coeff{i}{N+1}(:,2:end)).*Phase_part(:,2:end)),2);
        end
        n_order_length_vec = [n_order_length_vec,n_end-1];
        N = N +1;
%         clc; 
%         display(['Computing expansion coefficients for angular index: ',num2str(N)]);
        n = n_end;
    else
        % The case where the initial n isn't large enough -> Try again.        
        n = 2*n;
    end
end

vol_hat = zeros(gridSize,gridSize,gridSize,nVols);
vol_tmp = zeros(gridSize,gridSize,gridSize);
for i=1:nVols    
    vol_tmp(ball) = vol_mat_hat(:,i);
    vol_hat(:,:,:,i) = vol_tmp;
end

end

