function [coeff] = get_Gamma_Coeff(l,s,m,k,gamma)
%   Obtain a single gamma coefficient.
%   Boris Landa, 03.08.2018

coeff_l = gamma.coeff{m+1}(gamma.coeff_l_indices{m+1}==l,k+1);
coeff = coeff_l(s+1);

end

