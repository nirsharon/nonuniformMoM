function [val, grad] = LS_cost_w_grad(vec_AB, gamma, C_tensor, Gamma_mat , ...
                                        sign_mat, mu1_hat, mu2_hat, weight)

% A cost function for (complex) least squares optimization, with gradient!


% balancer between the two terms of mu1 and mu2
if isempty(weight)
    weight = 1;
end

% back to cell-array form
[A,B] = VecAB_to_A_B(vec_AB, gamma);

% initialization
K = size(gamma.coeff{1},2);

% first moment
m1     = FirstMoment_PSWF_v2(A, B, Gamma_mat, sign_mat);
inner1 = m1 - mu1_hat;
term1  = norm(inner1(:))^2;   % "vector" norm

% second moment
m2     = SecondMoment_PSWF_v2(A, B, gamma, C_tensor, Gamma_mat);
inner2 = m2 - mu2_hat;
term2  = norm(inner2(:))^2;   % "vector" norm

% function evaluation
val = term1 + weight*term2;  % just second moment at this point

% gradient evaluation, if needed
if nargout>1
    grad = get_gradient_PSWF(vec_AB, gamma, C_tensor, Gamma_mat, sign_mat,...
                                inner1, inner2, weight);
                            
    % filtering the (gradient of the) zero angular coefficient to be REAL
    % IF RUNNING GRAD CHECK MUST COMMENT IT!
    start_ind = 1;                        
    for j=1:length(A{1})
        s = size(A{1}{j},1);
        l_ind = start_ind:(start_ind+s-1);
        grad(l_ind) = real(grad(l_ind));
        start_ind = start_ind + numel(A{1}{j});
    end
end

end


