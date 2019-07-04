function [gamma_t,A_t] = gamma_truncate(gamma,A)
% Truncate the index sets associated with gamma and A to a maximal band L
% and a maximal angular frequency T, such that K(l)>=2l+1 and T(q)>=2l+1
% for all 0<=l<=L and 0<=|q|<=Q.
% Assumption:  T(l)>=K(l) for 0<=l<=L.

gamma_t=gamma;
A_t=cell(1);
for l=0:max(gamma.band_idx_3d)
    cutOff = nnz(gamma.band_idx_3d==l)<2*l+1;    
    if cutOff
        l_max = l-1;
        break;
    end
    A_t{1}{l+1} = A{1}{l+1};
end
gamma_t.band_idx_3d(gamma_t.band_idx_3d>l_max) = [];
gamma_t.ang_idx_2d(gamma_t.ang_idx_2d>l_max) = [];

for m=0:max(gamma.ang_idx_2d)        
    gamma_t.coeff{m+1}(nnz(gamma.coeff_l_indices{m+1}<=l_max)+1:end,:) = [];
    gamma_t.coeff_l_indices{m+1}(nnz(gamma.coeff_l_indices{m+1}<=l_max)+1:end,:) = [];             
end


end

