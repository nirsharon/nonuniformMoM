function [A_proj, B_proj] = Project_coefs_PSWF(A, B, SO_Array)
% This function project a given couple of cell arrays to a valid volume and
% distribution. The constraints are:
% 1 - Self refletion symmetry of A (representing a Fourier of real volume)
% 2 - Self refletion symmetry of B (representing a real distribution)
% 3 - Positivity of B over a given grid (relaxation of positivity)
% 4 - A unit inegration of B (B{1} = 1)

[AI,AE,Acon] = linear_cons_B2(numel(B)-1,SO_Array);
[B_proj, ~] = project_B_to_positive2(AI,AE,Acon,B);
factor1 = 1/B_proj{1};
for j=1:(length(B))
    B_proj{j} = factor1*B_proj{j};
end

% The constraint on A: A_l^(-m) = (-1)^{m+l}*conj(A_l^m). Namely,
% 1 - real(A_l^(-m)) = (-1)^{m+l}*real(A_l^m)
% =>  Aproj_l^m = 1/2( real(A_l^(-m)) + (-1)^{m+l}*real(A_l^m) )
% 2 - imag(A_l^(-m)) = (-1)^{m+l+1}*image(A_l^m)
% =>  Aproj_l^m = 1/2( imag(A_l^(-m)) + (-1)^{m+l}*imag(A_l^m) )
% and the negative m-s will follow

L = length(A);
A_proj = cell(L,1);
if L>0
    K = size(A{1},1);
    for l=0:(L-1)
        A_proj{l+1} = zeros(K,2*l+1);
        for m=0:l
            if m==0
                A_proj{l+1}(:,m+l+1) = ...
                    .5*((A{l+1}(:,l+1)) + (-1)^(l)*conj(A{l+1}(:,l+1)));
            else
                real_mean = .5*(real(A{l+1}(:,-m+l+1)) + (-1)^(m+l)*real(A{l+1}(:,m+l+1)));
                comp_mean = .5i*(imag(A{l+1}(:,-m+l+1)) + (-1)^(m+l+1)*imag(A{l+1}(:,m+l+1)));
                A_proj{l+1}(:,-m+l+1) = real_mean + comp_mean;
                A_proj{l+1}(:,m+l+1)  = (-1)^(m+l)*real_mean + (-1)^(m+l+1)*comp_mean;
            end
        end
    end
end


end

