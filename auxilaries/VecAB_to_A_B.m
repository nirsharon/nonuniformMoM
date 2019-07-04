function [A,B] = VecAB_to_A_B(vecAB, gamma)
%
% A utilizer -- parsing back the vector into the original cell arrays
%
% NS, Feb 19

L =  max(gamma.band_idx_3d)+1;

A_inner  = cell(1,L);
ind_list = 0;
for l=1:L
    current_s = nnz(gamma.band_idx_3d == (l-1));
    ind_list  = (ind_list(end)+1):(ind_list(end) + current_s*l );  
    A_inner{l} = reshape( vecAB(ind_list), current_s, l);
end

current = ind_list(end);
remain = length(vecAB) - current;
if remain>0
    Ppoly = [4,0,-1,-3*(remain)]; 
    P = roots(Ppoly);
    P = round(P(1));
    B = cell(P,1);
    B{1} = 1;
    starting_ind = current + 1;
    for p=2:P
        endin = (2*p-1)^2-1;
        B{p} = reshape(vecAB(starting_ind:(starting_ind+endin)),2*p-1,2*p-1);
        starting_ind = starting_ind + endin + 1;
    end
else
    B = 1;
end

A = {A_inner};
end


