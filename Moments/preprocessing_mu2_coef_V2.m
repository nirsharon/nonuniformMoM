function [C_tensor] = preprocessing_mu2_coef_V2(gamma, P, M)
%
% This function prepare the tensor of coefficients for mu2 and the lookup
% table of products of CG coefficients

% Input:
%   P = size(B,1), expansion length of the distribution
%   M = max(gamma.ang_idx_2d)+1 angular size of the moment
%
% NS, April 19

% the overall degree. Common to both 2D and 3D
L = max(gamma.band_idx_3d);

% just in cases where degree <= 30 NEED  TO UPDATE
assert(L<=30,'Degree higher than 30, check C-G table availablity');

% load the table of CG coefs
C = load('CGC_EZ_table_30'); C = C.C;
%C_array = zeros(M,M,L+1,L+1,L+1,L+1,P);

% index in vectorial B
ind_Bvec = @(p,u,v) (2*p-1)*(2*p)*(2*p+1)/6 + u + (v-1)*(2*p+1);
B_vec_length   = (2*P-1)*(2*P)*(2*P+1)/6;

% initialization
C_size   = (L+1)^2; % old: (L^2+3*L+2)/2;
C_tensor = cell(M, M);

% index in C
l_max_ind = @(l) (l+1)^2; %  old: (l^2+3*l+2)/2;

% main loop
for m1 = 0:M-1
    for m2 = 0:M-1
        % initialize current matix of coefficients
        C_tensor{m1+1,m2+1} = zeros(C_size^2,B_vec_length);
        
        for l1 = m1:L
            l1_ind = (l_max_ind(l1-1) +1):l_max_ind(l1);

            for l2 = m2:L
                l2_ind = (l_max_ind(l2-1) +1):l_max_ind(l2);

                % the row indices for the tensor
                combined_l_ind = C_size*(l2_ind-1) + l1_ind.';
                
                % inner loop
                for p = abs(l1-l2):min(l1+l2,(P-1))
                    if (abs(m1-m2)>p)
                        continue
                    else
                        % get the relevant n1 and n2 indices
                        ind = (-l2:l2)-(-l1:l1)';
                        lin_ind_ar = 1:numel(ind);
                        nnz_ind = lin_ind_ar(abs(ind)<=p);
                        [ind_n1,ind_n2] = ind2sub(size(ind),nnz_ind);
                        n1_arr = -l1:l1;
                        n2_arr = (-l2:l2);
                        
                        % get the relevant n1 and n2 indices
                        mn1 = max(ind_n1);
                        mn2 = max(ind_n2);
                        c_n1_rang = min(n1_arr(ind_n1)):max(n1_arr(ind_n1));
                        c_n2_rang = min(n2_arr(ind_n2)):max(n2_arr(ind_n2));
                        
                        C_here = C{l1+1,p+1,l2+1}(c_n1_rang + l1+1, c_n2_rang + l2+1)...
                            *C{l1+1,p+1,l2+1}(m1+1+l1,m2+1+l2)/(2*l2+1);
                        
                        C_array_p = reshape(C_here, mn1, mn2);
                        
                        %C_array_p = reshape(C_array(l1+1, l2+1, m1+1, m2+1 ,1:mn1, 1:mn2, p+1), mn1, mn2);
                        C_array_p(lin_ind_ar(abs(ind)>p)) = 0;   % for double check..
          
                        % getting the right B indices, namely columns
                        B_ind = ind_Bvec(p, n2_arr(ind_n2)-n1_arr(ind_n1) +p+1, (m2-m1) +p+1);
                        v     = C_array_p(nnz_ind);
                        l_ind = combined_l_ind(nnz_ind);
                        
                        % fill the matrix
                        linearInd = sub2ind(size(C_tensor{m1+1,m2+1}), l_ind(:), B_ind(:));
                        C_tensor{m1+1,m2+1}(linearInd) = v(:);                       
                        
%                         
%                         % the coefficients in the C saved array
%                         % the required product
%                         C_array(l1+1, l2+1, m1+1, m2+1 ,c_n1_rang + l1+1, c_n2_rang + l2+1, p+1) = ...
%                             C{l1+1,p+1,l2+1}(c_n1_rang + l1+1, c_n2_rang + l2+1)...
%                             *C{l1+1,p+1,l2+1}(m1+1+l1,m2+1+l2)/(2*l2+1);
                    end
                end  % p loop
                
            end      % l2 loop
        end          % l1 loop
    end              % m2 loop
end                  % m1 loop


end

