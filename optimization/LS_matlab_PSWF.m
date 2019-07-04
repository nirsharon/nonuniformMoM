function [A_est, B_est] = LS_matlab_PSWF(initial_guess, m1, m2, gamma, ...
    C_tensor, Gamma_mat , sign_mat, iters)

options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true,'Display','iter');
options.OptimalityTolerance = 1e-16;
options.StepTolerance       = 1e-16;
if  ~exist('iters','var')
    options.MaxIterations       = 200;
else
    options.MaxIterations       = iters;
end

options.SubproblemAlgorithm = 'factorization';

x0    = [real(initial_guess); imag(initial_guess)];  
weight = nnz(m1)/nnz(m2)*50;
fobj = @(x) get_real_costgrad(gamma, C_tensor, Gamma_mat , ...
                                        sign_mat, m1, m2, weight, x);
% run
x    = fminunc(fobj, x0, options);



% summary
N     = numel(x)/2;
vecAB = x(1:N) + 1i*x(N+1:2*N);
[A_est, B_est] = VecAB_to_A_B(vecAB, gamma);

end


% from complex to Matlab format
function [cost, grad] = get_real_costgrad(gamma, C_tensor, Gamma_mat , ...
                                        sign_mat, m1, m2, weight, x)

N  = numel(x)/2;

% translate to complex variable
x_complex = x(1:N)+1i*x(N+1:2*N);

% apply the function
switch nargout
    case 1
        [cost] = LS_cost_w_grad(x_complex, gamma, C_tensor, Gamma_mat , ...
                                        sign_mat, m1, m2, weight);
    case 2
        [cost, grad_complex] = LS_cost_w_grad(x_complex, gamma, ...
            C_tensor, Gamma_mat, sign_mat, m1, m2, weight);
        % return the grad
        grad = [real(grad_complex); imag(grad_complex)];
end



end