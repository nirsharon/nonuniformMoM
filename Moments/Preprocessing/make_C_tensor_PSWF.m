function [C_tensor] = make_C_tensor_PSWF(gamma, C_array, P)
%
%  This function calculates the tensor coefficient for mu2, prolates
%  version
%
%           Fixed version, April 19
%
% Input:
%    gamma   -- a struct that includes all info regarding the conversion between
%                PSWF coefficients of 2D and 3D objects
%    C_array -- array of clebsch-gordon coefficients (preprocessing)
%    P = size(B,1), expansion length of the distribution
%
% Output:
%    C_tensor -- the tensor, st for g = Gamma*A, the (m1,m2) entry of mu2 is
%         g*reshape(C_tensor*vec_B, C_size, C_size)*g'
%
% NS, Mar 19
%------------------------------------------
% NOTE: apply "FUTURE CHANGE" when C array of CG is fully avaiable                      
%------------------------------------------

% FUTURE CHANGE (uncomment next two lines)
% load and initilize
% C = load('CGC_EZ_table_30'); C = C.C;

% the overall degree. Common to both 2D and 3D
L = max(gamma.band_idx_3d);

% total sizes of the moment
sizeM = max(gamma.ang_idx_2d)+1 ;
C_size = (L+1)^2; % old: (L^2+3*L+2)/2;

% index in vectorial B
ind_Bvec = @(p,u,v) (2*p-1)*(2*p)*(2*p+1)/6 + u + (v-1)*(2*p+1);
B_vec_length   = (2*P-1)*(2*P)*(2*P+1)/6;

% initialization
C_tensor = cell(sizeM,sizeM);

% index in C
l_max_ind = @(l) (l+1)^2; %  old: (l^2+3*l+2)/2;

% main loop
for m1 = 0:sizeM-1
    for m2 = 0:sizeM-1
        
        % initialize current matix of coefficients
        C_tensor{m1+1,m2+1} = zeros(C_size^2,B_vec_length);
        
        % run over l-s
        for l1 = m1:L
            %l1_ind = (((l1-1)^2+3*(l1-1)+2)/2+1):((l1^2+3*l1+2)/2);
            l1_ind = (l_max_ind(l1-1) +1):l_max_ind(l1);
            
            for l2 = m2:L
                %l2_ind = (((l2-1)^2+3*(l2-1)+2)/2+1):((l2^2+3*l2+2)/2);
                l2_ind = (l_max_ind(l2-1) +1):l_max_ind(l2);

                % the row indices
                combined_l_ind = C_size*(l2_ind-1) + l1_ind.';
                
                % inner loop
                for p=abs(l1-l2):min(l1+l2,(P-1))
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
                        mn1 = max(ind_n1);
                        mn2 = max(ind_n2);
                        
                        % the relevant coefficients
                        
                        % FUTURE CHANGE (uncomment next two lines)
                        %c_prod = C{l1+1,p+1,l2+1}(c_ind_n1,c_ind_n2)*C{l1+1,p+1,l2+1}(m1+1+l1,m2+1+l2)/(2*l2+1);  
                        %C_array_p = reshape(c_prod, mn1, mn2);

                        c_n1_rang = min(n1_arr(ind_n1)):max(n1_arr(ind_n1));
                        c_n2_rang = min(n2_arr(ind_n2)):max(n2_arr(ind_n2));
                        C_array_p = reshape(C_array(l1+1, l2+1, m1+1, m2+1 ,c_n1_rang + l1+1, c_n2_rang + l2+1, p+1), mn1, mn2);
                        
                        %C_array_p = reshape(C_array(l1+1, l2+1, m1+1, m2+1 ,1:mn1, 1:mn2, p+1), mn1, mn2);
                        C_array_p(lin_ind_ar(abs(ind)>p)) = 0;   % for double check..
                        
                        % getting the right B indices, namely columns
                        B_ind = ind_Bvec(p, n2_arr(ind_n2)-n1_arr(ind_n1) +p+1, (m2-m1) +p+1);
                        v     = C_array_p(nnz_ind);
                        l_ind = combined_l_ind(nnz_ind);
                        
                        % fill the matrix
                        linearInd = sub2ind(size(C_tensor{m1+1,m2+1}), l_ind(:), B_ind(:));
                        C_tensor{m1+1,m2+1}(linearInd) = v(:);
                    end
                end       % inner loop
            end           % l2 loop
        end               % l1 loop
    end                   % m2 loop
end                       % m1 loop

end


