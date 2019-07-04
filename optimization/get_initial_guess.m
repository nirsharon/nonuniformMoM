function [x0] = get_initial_guess(P, AI, AE, Acon, M)
% Input:
%   P - length of distribution expansion
%   AI, AE, Acon - the constraints on the distribution
%   M - the required length, that is length(vec_AB_GT)
% Output:
%   x0 - random initial guess. The second part has the structure of
%   a distribution
%
% NS, May 19

% distributon expansion coefficients
B_init = cell(P,1);
pert_mag = .1;
B_init{1} = 1;
for j = 2:P
    B_init{j} = pert_mag*(rand(2*j-1)*.1 + rand(2*j-1 )*.1*1i)/(2*j-1) ;
end
[B_init, ~] = project_B_to_positive2(AI, AE, Acon, B_init);
scl   = 1/B_init{1};
for i = 1:numel(B_init)
    B_init{i} =  scl*B_init{i};
end

% construct the initial guess all together
x0 = rand(M,1);
vb = FromCellArr2Vec({1},B_init);
vb = vb(2:end);
x0((end-length(vb)+1):end) = vb;

end

