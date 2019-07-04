function [A_coef_mu2] = get_A_coefs_mu2_V2(A, B, Gamma_mat, C_tensor)

% The coefficients of A in the moment
%
% NS, April 19

% total sizes 
sizeM  = length(Gamma_mat);
sizeK  = size(Gamma_mat{1},1);
L      = sizeM-1;
C_size = (L+1)^2;

% vectorized versions
vec_A   = A_B_to_VecAB(A, [], sizeM-1);
vec_B   = FromCellArr2Vec({1},B);
lengthA = length(vec_A);

% index matrix
ind_mat_real = [];
ind_mat_conj = [];
A_mat        = [];

% prepare index matrices
startA = 0;
for m = 0:sizeM-1
    Al_s_n      = A{1}{m+1}(:,1:(m+1));
    current_blk = [ conj(Al_s_n(:,end:-1:2))*diag((-1).^(-m:1:-1).') , Al_s_n];
    A_mat = blkdiag(A_mat,current_blk);
    
    v = (1:numel(Al_s_n)) + startA;
    startA = v(end);
    
    current_Aind = reshape(v,size(Al_s_n));
    blk_ind_real = [ zeros(size(current_Aind(:,end:-1:2))) , current_Aind];
    ind_mat_real = blkdiag(ind_mat_real, blk_ind_real);
    blk_ind_conj = [ current_Aind(:,end:-1:2)*diag((-1).^(-m:1:-1).') ...
        , zeros(size(current_Aind))];
    ind_mat_conj = blkdiag(ind_mat_conj, blk_ind_conj);
end

% initializing
A_coef_mu2      = zeros(sizeM, sizeK, sizeM, sizeK, lengthA);

% main loop
for j=1:lengthA
    [x,y] = ind2sub(size(ind_mat_real),find(ind_mat_real==j));
    [a,b] = ind2sub(size(ind_mat_conj),find(abs(ind_mat_conj)==j));
    signj = sign(ind_mat_conj(a,b));
    
    for m1 = 0:sizeM-1
        A1 = Gamma_mat{m1+1}*A_mat;
        K1 = size(Gamma_mat{m1+1},1);
        for m2 = 0:sizeM-1
            A2 = Gamma_mat{m2+1}*A_mat;
            K2 = size(Gamma_mat{m2+1},1);
            
            C = reshape(C_tensor{m1+1,m2+1}*vec_B, C_size, C_size);
            
            part1 = Gamma_mat{m1+1}(:,x)*C(y,:)*A2';
            part2 = A1*C(:,y)*(Gamma_mat{m2+1}(:,x)');
            current_block = part1 + part2;
            
            if ~isempty(signj)
                part1_conj = signj*Gamma_mat{m1+1}(:,a)*C(b,:)*A2';
                part2_conj = signj*A1*C(:,b)*(Gamma_mat{m2+1}(:,a)');
                current_block = current_block + conj(part1_conj+part2_conj);
            end
            
            A_coef_mu2(m1+1, 1:K1, m2+1, 1:K2, j) = current_block;
            
        end
    end
    
end

end


