function [coeffs] = pswf_t_f_3d(vol_arr, beta, delta)
%   This is the forward 3D PSWF transform - from volumes sampled on the Cartesian
%   grid to the 3D PSWF expansion coefficients.
%   The 3D PSWFs are given by \psi_{N,n,m}(r,\theta,\phi) = R_{N,n}(r) * Y_N^m(\theta,\phi),
%   where R_{N,n}(r) are radial functions and Y_N^m(\theta,\phi) are the spherical harmonics.
%   Input:  vol   :     Volume functions (real-valued) to be expanded. The fourth dimension enumerates over different volumes, where the
%                       first three dimensions correspond to an equally-spacsed grid of points.
%           beta:       Bandlimit ratio relative to the Nyquist rate, between 0 and 1.
%           delta:      Truncation parameter, between 0 and 1, which controls the length of the
%                       expansion and the approximation error. Smaller values (closer to zero) guarantee smaller errors, yet longer expansions, and vice-versa. 
%                       Note: Due to numerical considerations, do not exceed [1E-6,1-1E-6].
%   Output: coeffs:     PSWF expansion coefficients, provided as a cell
%                       array enumrated by the index of the volume function. Within each cell value, there is another cell array enumerated by the angular index N (equall to \ell in the classical spherical harmonics notation), where each cell contains a
%                       matrix of expansion coefficients such that the first dimension (rows)
%                       correspond to the radial index n, and the second dimension (columns)
%                       correspond to the spherical harmonic index m (which is between 0 and N since the volume is assumed to be real-valued).
%
%
%   Boris Landa, 04.12.2018
%
%% Compute expansion coeffs for every angular frequency
[sizeX,~,~,nVols] = size(vol_arr);
L = floor(sizeX/2);
if mod(sizeX,2)==0
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

% nVols = size(vol_arr,4);
vol_mat = zeros(nnz(ball),nVols);
for i=1:nVols
    currVol = vol_arr(:,:,:,i);
    vol_mat(:,i) = currVol(ball);
end

c = beta * pi*L;

% coeffs = {};
coeffs = cell(1,nVols);
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
        PSWF_R = PSWF_R*diag(abs(alpha).^2)*(beta/2)^3;
        PSWF_R = PSWF_R(:,1:n_end-1);
        Phase_part = mySph_v2(N,0:N,theta_b,phi_b);     % Complex valued basis functions
        for i=1:nVols
            vol_phase = bsxfun(@times,conj(Phase_part),vol_mat(:,i));
            coeff = PSWF_R' * vol_phase;
            for m=0:N            
                if (m==0)
                    coeffs{i}{N+1} = coeff(:,m+1);
                else
                    coeffs{i}{N+1} = [coeffs{i}{N+1} coeff(:,m+1)];
                end
            end        
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

end

