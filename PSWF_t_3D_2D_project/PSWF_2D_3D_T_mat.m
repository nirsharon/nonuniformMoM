function [gamma] = PSWF_2D_3D_T_mat(c, delta, eps_p)
%   Compute mapping coefficients between volume 3D PSWF coefficients and 2D
%   projection image PSWF coefficients.
%   Input:  c:          Bandlimit in [rad/sec], i.e. c = pi*f_s for Nyquist, where f_s is the sampling frequency in [Hz]. 
%           delta:      Truncation parameter, between 0 and 1, which controls the length of the
%                       expansion in 3D spatial domain. Values closer to 1 lead to shorter expansions aimed for localized volume functions. 
%           eps_p:      Prescribed accuracy (in squared error norm) for 2D projection image approximation by 2D PSWFs. This controlls the expansion length in 2D. 
%   Output: gamma:      Gamma coefficients arranged into a cell array of
%   matrices. Every cell corresponds to a different angular index m in the
%   2D PSWf expansion. The rows of a matrix in each cell correspond
%   to different pairs (l,s), where l are band indices (running from m to the maximal band index) and s are radial indices,
%   and the columns correspond to different radial indices k of the 2D PSWF
%   expansion.
%
%   Boris Landa, 03.08.2018
%
%% Compute 3D PSWFs radial part
[r,w]=lgwt(1e3,0,1);
N = 0;
n = ceil(2*c/pi);
F_tot = [];
band_idx = [];
alpha = [];
while(1)
    [F,alpha_ls] = PSWF_ND(N,n,c,1,r,eps);
    lambda = (c/(2*pi))^3 * abs(alpha_ls).^2;
    n_end = find(lambda<delta,1,'first');
    
    if (n_end < 2)
        break;
    end
        
    if (~isempty(n_end))
%         Phase_part = mySph_v2(N,0:N,theta_b,phi_b);     % Complex valued basis functions                   
        F_tot = [F_tot F(:,1:n_end-1)];
        alpha = [alpha alpha_ls(1:n_end-1).'];
        band_idx = [band_idx repmat(N,1,n_end-1)];
        N = N +1;
        n = n_end;
    else
        % The case where the initial n isn't large enough -> Try again.        
        n = 2*n;
    end
end

%% Compute 2D PSWFs radial part and mapping coefficients
delta_prime = 1e-6;
N = 0;
n = ceil(2*c/pi);
f_tot = [];
ang_idx = [];
beta = [];
coeff ={};
coeff_indices = {};
truncation_tail = [];
ang_idx_2d = [];
while(1)
    [f,beta_mk] = PSWF_ND(N,n,c,0,r,eps);
    mu = (c/(2*pi))^2 * abs(beta_mk).^2;
    n_end = find(mu<delta_prime,1,'first');
    
    if (n_end < 2)||(N>max(band_idx))
        break;
    end
        
    if (~isempty(n_end))
%         Phase_part = mySph_v2(N,0:N,theta_b,phi_b);     % Complex valued basis functions                   
        f_tot = [f_tot f];
        beta = [beta beta_mk.'];
        ang_idx = [ang_idx repmat(N,1,size(f,2))];
        % - Compute coefficients matrix         
        Ysh=[];
        for l=N:max(band_idx)
            Ysh = [Ysh; repmat(mySph_v2(l,N,pi/2,0),nnz(band_idx==l),1)];
        end
        H = F_tot(:,band_idx>=N)' * diag(w.*r) * f_tot(:,ang_idx==N);
%         gamma_all = 2*pi/c*((Ysh./(alpha(band_idx>=N).'))*beta_mk.').*H;
        gamma_all = (2*pi)^(3/2)/c*((Ysh./(alpha(band_idx>=N).'))*beta_mk.').*H;
        truncation_tail(N+1,:) = sum(flipud(cumsum(flipud(repmat(2*band_idx(band_idx>=N).'+1,1,size(gamma_all,2)).'.*abs(gamma_all).^2.'))).',1);
        coeff{N+1} = gamma_all(:,truncation_tail(N+1,:)>eps_p);
        coeff_indices{N+1} = band_idx(band_idx>=N);        
        ang_idx_2d = [ang_idx_2d repmat(N,1,size(coeff{N+1},2))];
%         truncation_tail = fliplr(cumsum(fliplr(max(abs(gamma)).^2 * sum(2*band_idx(band_idx>=N).'+1))));
        N = N +1;
%         n = n_end;        
    else
        % The case where the initial n isn't large enough -> Try again.        
        n = 2*n;
    end
end

gamma.band_idx_3d = band_idx;
gamma.ang_idx_2d = ang_idx_2d;
gamma.coeff = coeff;
gamma.coeff_l_indices = coeff_indices;

%% Display residue error for truncation
% figure; 
% semilogy([truncation_tail.']); grid on; axis([-inf 25 1e-16 inf]); 
% ylabel('$\epsilon_m(k)$','Interpreter','Latex');
% xlabel('$k$','Interpreter','Latex')

end

