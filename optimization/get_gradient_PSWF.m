function [grad] = get_gradient_PSWF(x, gamma, C_tensor, Gamma_mat, sign_mat,...
                                     remainder1, remainder2, weight)
% 
% This function calculates the gradient of the least squares cost function;
% 
% Input
%     x -- the vector of unknowns (A plus B))
% gamma -- the structure of coordinate changing.
% C_tensor,Gamma_mat, sign_mat -- the auxilary matrices for the moments
%                                 calculations
% remainder1, remainder2       -- the inner parts of the terms of the first 
%                                 and second moments, respectively.
% weight -- the balancer between the two terms

%
% NS, April 19 (fixed)

% back to cell arrays
[A, B] = VecAB_to_A_B(x, gamma);

% max degree, total size
max_deg = min(length(A{1})-1, size(B,1)-1);  
sizeM   = max(gamma.ang_idx_2d)+1 ;
L       = max(gamma.band_idx_3d);
C_size  = (L+1)^2; 

% toward generating index matrix
A_ind = A;
all_A_dependent_ind = 0;
just_positive_A_ind = 0;
startA = 0;
for m = 0:sizeM-1
    v = (1:numel(A_ind{1}{m+1})) + startA;
    A_ind{1}{m+1} = reshape(v,size(A_ind{1}{m+1}));
    startA = v(end);
    all_A_dependent_ind = all_A_dependent_ind + ...
        numel([A_ind{1}{m+1},A{1}{m+1}(:,end:-1:2)]);
    just_positive_A_ind = just_positive_A_ind + numel(A_ind{1}{m+1});
end

% prepare A to blkdiag
Mat = zeros(all_A_dependent_ind,just_positive_A_ind);
% move A to blkdiag and save the indix positions
A_mat       = [];
mat_ind_all = [];
mat_ind     = [];
current_ind = 1;
for m = 0:sizeM-1
    Al_s_n      = A{1}{m+1}(:,1:(m+1));
    current_blk = [ conj(Al_s_n(:,end:-1:2))*diag((-1).^(-m:1:-1).') , Al_s_n];
    current_Aind = A_ind{1}{m+1}(:,1:(m+1));
    blk_ind = [ current_Aind(:,end:-1:2)*diag((-1).^(-m:1:-1).') , current_Aind];
            
    A_mat = blkdiag(A_mat,current_blk);

    % prepare a matix that will save the sign flipping
    current_ind = current_ind:(current_ind + numel(current_blk)-1);
    lin_ind_mat = sub2ind(size(Mat), current_ind(:) ,abs(blk_ind(:)));
    Mat(lin_ind_mat) = sign(blk_ind(:));
    current_ind = current_ind(end) + 1;
    
    % matrices that will save the relevant indices to extract
    relevant_ind = zeros(size(current_blk));
    relevant_ind(:,(end-size(Al_s_n,2)+1):end) = 1;
    mat_ind     = blkdiag(mat_ind,relevant_ind);
    mat_ind_all = blkdiag(mat_ind_all,ones(size(current_blk)));   
end
diff_ind = mat_ind_all(:)-mat_ind(:);
conj_ind = (diff_ind(mat_ind_all(:)>0));

% vectorize B
vec_B = FromCellArr2Vec({1},B);

% main loop
mu1_grad = 0;
part1 = 0; 
part2 = 0;
for m1 = 0:sizeM-1
    A1 = Gamma_mat{m1+1}*A_mat;
    K1 = nnz(gamma.ang_idx_2d==m1);
    % first moment part
    if m1<=max_deg
        B_part = 2*(Gamma_mat{m1+1}*A_mat*sign_mat{m1+1})'*remainder1(m1+1,1:K1).'; 
        A_part = 2*Gamma_mat{m1+1}'*remainder1(m1+1,1:K1).'*(sign_mat{m1+1}*vec_B)';        
        
        full_vec_coef = (A_part(mat_ind_all(:)>0));
        conj_part = conj(Mat.'*(full_vec_coef.*conj_ind));

        A_all = A_part(mat_ind(:)>0) +  conj_part;
        
       % mu1_grad = mu1_grad + [A_part(mat_ind(:)>0); B_part(2:end)];        
        mu1_grad = mu1_grad + [A_all; B_part(2:end)];        
    end
    
    for m2 = 0:sizeM-1
        
        A2 = Gamma_mat{m2+1}*A_mat;
        K2 = nnz(gamma.ang_idx_2d==m2);

        % first part - the volume
        r = squeeze(remainder2(m1+1,1:K1, m2+1,1:K2));
        Big_C = reshape(C_tensor{m1+1,m2+1}*vec_B, C_size, C_size);
        md    = 2*Gamma_mat{m2+1}'*r'*A1*Big_C + 2*Gamma_mat{m1+1}'*r*A2*Big_C';        
        
        full_vec_coef = (md(mat_ind_all(:)>0));
        conj_part = conj(Mat.'*(full_vec_coef.*conj_ind));
        
        A_all = md(mat_ind(:)>0) +  conj_part;
        
        % summary of A part
        part1 = part1 + A_all;
        
        % second part - the distribution              
        B_part = 0;
        for k1 = 1:K1
            for k2 = 1:K2
                vec_a = kron(A1(k1,:)',A2(k2,:));  % can be done off line ?
                C_a   = ( vec_a(:)'*C_tensor{m1+1,m2+1} )';
                B_part = B_part + 2*(r(k1,k2))*C_a;
            end
        end
        part2  =  part2 + B_part(:);
                
    end
end

% conclusion
grad = mu1_grad + weight*[part1(:); part2(2:end)];

end

