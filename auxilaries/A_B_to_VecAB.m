function [vecAB] = A_B_to_VecAB(A, B, idx_3d)
%
% Converting the coefficient cells A and B to a vector
%
%  run as "A_B_to_VecAB(A, B, gamma.band_idx_3d)"
%
% NS, Feb 19

L = length(A{1});
P = size(B,1);

% A part (volume)
%vecAB = []; 
vecAB = zeros( sum(idx_3d)+ size(idx_3d,2) + (2*P-1)*(2*P+1)*P/3 -1 ,1);
ind_list = 0;
for l=1:L
    current_v = A{1}{l}(:);
    ind_list  = (ind_list(end)+1):(ind_list(end) + numel(current_v) );
    vecAB(ind_list) = current_v;
end

% B part (distribution)
current_ind    = ind_list(end)+1; 
b_range        = current_ind:(current_ind+(2*P-1)*(2*P+1)*P/3 - 2);
vecAB(b_range) = 0;
for p=2:P
    endin = (2*p-1)^2-1;
    vecAB(current_ind:(current_ind+endin)) = B{p}(:);
    current_ind = current_ind + endin + 1;
end

end

