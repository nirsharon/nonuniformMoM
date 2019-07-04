function [R,alpha_N] = PSWF_3D(N,n,c,radial_eval_pts,approx_err)
%% Generate points for numerical integration
[x,w]=lgwt(20*max(N,n),0,1);    % Legendre-Gauss nodes for numerical integration   
x = flipud(x);
w = flipud(w);

p = 1;    % dimension - D = 2 + p 
N_tilde = N + p/2;

%% Calculate approx_len
% - Definition for the approximation function of the d's decay
d_decay_approx_fun = @(N,n,c,j) c^2./(16*(j.^2 + j*(2*n+N+1)) - c^2);

% - Finding the first point of the decay according to the approximation
first_idx_for_decrease = ceil((((2*n+N+1)^2+c^2/2)^0.5 - (2*n+N+1))/2);     % Analytic solution for d_decay_approx_fun<1 for every j>first_idx_for_decrease

% - Finding the number of points needed to achieve the chosen error
d_approx = d_decay_approx_fun(N_tilde,n,c,first_idx_for_decrease);
d_decay_idx_counter = first_idx_for_decrease;

while(d_approx>approx_err)
    d_decay_idx_counter = d_decay_idx_counter + 1;
    d_approx = d_approx * d_decay_approx_fun(N_tilde,n,c,d_decay_idx_counter);
end

approx_len = n + 1 + d_decay_idx_counter;

%% Funcion definitions
Pn = @(n,alpha,beta,x) j_polynomial(length(x),n,alpha,beta,x);

%% Generate Matrix for calculation of d_k coefficients
B_N = Generate_BN_mat(N_tilde,c,approx_len);
[d_vec, xi] = eig(B_N);

%% Sort eigenvalues by descending order
[xiVec,idx] = sort(diag(xi),'descend');
d_vec = d_vec(:,idx);
for i=1:approx_len
    xi(i,i) = xiVec(i);
end

%% Compare true decay of d's with approximation function
% d_decay_compare = [abs(d_vec((n+2):end,n+1)./d_vec((n+1):(end-1),n+1)),d_decay_approx_fun(N,n,c,0:(approx_len-2-n)).'];
% d_approx = cumprod([ones(first_idx_for_decrease,1);d_decay_approx_fun(N,n,c,first_idx_for_decrease:(approx_len-1-n)).']);
% d_decay_compare = [abs(d_vec((n+1):end,n+1)),d_approx];
% plot(d_decay_compare); grid on; title(['N=',num2str(N),', n=',num2str(n),', c=',num2str(c)] );

%% Generate T_n basis functions on grid points ,and the derivatives of T_n.
j=0:(approx_len-1);

% - Calculation of functions
T_x_mat = @(x,N,j,approx_len) ((x.^(N+1/2))*((2*(2*j+N+1)).^(1/2))) .* Pn(approx_len-1,N,0,1-2*x.^2);
% - Calculation of derivatives
T_x_derivative_mat = ((x.^(N_tilde+3/2))*(-2*(N_tilde+j+1).*((2*(2*j+N_tilde+1)).^(1/2)))).*[zeros(length(x),1),Pn(approx_len-2,N_tilde+1,1,1-2*x.^2)] + (N_tilde+1/2)*((x.^(N_tilde-1/2))*((2*(2*j+N_tilde+1)).^(1/2))).*Pn(approx_len-1,N_tilde,0,1-2*x.^2);

%% Calculate phi - Generalized PSWF's ,and the derivatives of phi.
phi = T_x_mat(x,N_tilde,j,approx_len) * d_vec(:,1:(n+1));
% phi = phi .* (ones(length(x),1)*sign(phi(300,:)));    %Fix the signs so that all functions begin positive.

phi_derivatives = T_x_derivative_mat * d_vec(:,1:(n+1));
% phi_derivatives = phi_derivatives .* (ones(length(x),1)*sign(phi(300,:)));    %Fix the signs of the derivatives to be consistant with the signs of the functions.

%% Calculate original problem eigenvalues gamma_N,n

% - First (largest) eigenvalue numerical calculation
[max_phi_val,max_phi_x_idx] = max(abs(phi(:,1)));
x_for_calc = x(max_phi_x_idx);

K_operator = @(nu,x) besselj(nu,x).*sqrt(x);
Right_hand_side_integral = (w.*K_operator(N_tilde,c*x_for_calc*x)).' * phi(:,1);
gamma_N_1 = Right_hand_side_integral/phi(max_phi_x_idx,1);

% - Calculate the rest of the eigenvalues based on the first
upper_integral_values = diag( (((w.*x)*ones(1,n)).*phi_derivatives(:,1:(end-1))).' * phi(:,2:end) );
lower_integral_values = diag( (((w.*x)*ones(1,n)).*phi(:,1:(end-1))).' * phi_derivatives(:,2:end) );
gamma_N = [gamma_N_1 ; gamma_N_1 * cumprod(upper_integral_values./lower_integral_values)];

alpha_N = gamma_N*(2*pi)^(1+p/2)*((1i)^N) / (c)^((p+1)/2);

%% Test orthonormality of phi's
% Inner_prod_matrix = ((w*ones(1,n+1)).*phi).' * phi;
% figure; bar(abs(eig(Inner_prod_matrix)));

%% Verify numerically integral equation (test is effective only for eigenvalues in transitional region, otherwise tests only eigenvector subspace accuracy)
% numicIntEqErr = K_operator(N_tilde,c*x*x.')*diag(w)*phi - phi*diag(gamma_N);
% max(max(abs(numicIntEqErr)))

%% Compute PSWF's on the grid inside the disk
j=0:(approx_len-1);
T_radial_part_mat = @(x,N,j,approx_len) ((x.^(N-p/2))*((2*(2*j+N+1)).^(1/2))) .* Pn(approx_len-1,N,0,1-2*x.^2); % - These are the ordinairy T's divided by x^(1/2).
% d_vec = bsxfun(@times,d_vec,sign(sign(diag(d_vec).')+0.5));    %Fix the signs so that all functions begin positive.
maxVal = max(d_vec);
minVal = min(d_vec);
signVec = ones(1,size(d_vec,2));
signVec(abs(maxVal)<abs(minVal)) = -1;
d_vec = bsxfun(@times,d_vec,signVec);    % Fix the signs for consistency between runs (the eigenvectors can have both signs)
R = T_radial_part_mat(radial_eval_pts,N_tilde,j,approx_len) * d_vec(:,1:(n+1));
% R_radial_part_mat = R_radial_part_mat .* (ones(length(radial_eval_pts),1)*sign(phi(100,:)));    %Fix the signs so that all functions begin positive.
end