% script name: "test_Jacobian_coefs_V3_JK"
%
% In this version we construct the system with full vec of volume
% coefficients (no conjugated part is considered)

clear; clc;
fprintf('\n \n \n \n');

% generate GT and gamma 
% tic
% fprintf('Generating ground truth data..');
% P        = 3;
% S        = load('SO3_fifteen.mat');
% SO_grid  = S.SO3;
% gridSize = 91; % play with this...
% delta    = .9999999;     % .999999999999;
% [A, B, vol, gamma] = createGTdata(P, gridSize, SO_grid, delta);


gridSize = 21; % play with this...
radius = floor(gridSize/2);
beta=1;
c = beta*pi*radius;
delta    = .99;     % .999999999999;
eps_p = 1e-3;
gamma = PSWF_2D_3D_T_mat(c, delta, eps_p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Joe's try at the Jacobian.breakpoint/comment this on 1st run to get gamma
P=8;
L=4;

ind3D = gamma.band_idx_3d;
S = zeros(1,L+1);
for i=1:L+1 %really 0 to L
    S(i)=sum(ind3D(:)==i-1);
end
ind2D = gamma.ang_idx_2d;
T = zeros(1,L+1);
for i=1:L+1 %really 0 to L
    T(i)=sum(ind2D(:)==i-1);
end
Gamma=cell(L+1,L+1);%gamma^q,ell with q>=0.  Use gamma^q,ell = gamma^-q,ell..

%joe=gamma.band_idx_3d;
%L=10
%thirdD = zeros(1,L+1);
%for i=1:L+1
%thirdD(1,i)=sum(joe(:)==i-1);
%end



% P=3;%my P is different from above by 1 less, mine is actual max degree, uniform <-> P=0 
% L=4;
% S=[11 11 10 10 9];%radial bandlimits on volume, manual input at this point
% T=[13 13 12 12 11];%radial bandlimits on image, manual input, try T>=S now
% Gamma=cell(L+1,L+1);%gamma^q,ell with q>=0.  Use gamma^q,ell = gamma^-q,ell...
%                                   % maybe up to inconsequential sign (-1)^q

for q=1:L+1%really 0 to L
   for ell=q:2:L+1%only need to look at ell>=|q| and ell=q (mod2)
       Gamma{q,ell}=zeros(T(q),S(ell));
       for s=1:S(ell)
           for t=1:T(q)
               Gamma{q,ell}(t,s)=get_Gamma_Coeff(ell-1,s-1,q-1,t-1,gamma);
           end
       end
   end   
end%done extracting everything need (gamma's, really eta's) from Boris code



sizeA = 0;
for ell = 1:L+1 %really 0:L 
    sizeA = sizeA + (2*ell-1)*S(ell);
end


sizeB = 0;
for p = 0:P
    sizeB = sizeB + (2*p+1)^2;
end %count for totally non-uniform

sizeM1 = 0;
for q1=-min(L,P):min(L,P)
    sizeM1 = sizeM1 + T(abs(q1)+1);
end %count for totally non-uniform

sizeM2 = 0;
for q1 = -L:L
    for q2 = -L:q1
        if (abs(q1+q2)<=P)
            if q1 > q2 
                sizeM2 = sizeM2 + T(abs(q1)+1)*T(abs(q2)+1);
            else %q1=q2
                sizeM2 = sizeM2 + nchoosek(T(abs(q1)+1)+1,2);
            end
        end
    end
end %count for totally non-uniform

A = cell(L+1,1);
realAlpha = cell(L+1,1);
for ell=0:L
    A{ell+1}= zeros(S(ell+1),2*ell+1);
    realAlpha{ell+1} = 2*rand(S(ell+1),2*ell+1)-ones(S(ell+1),2*ell+1);%hopefully iid uniform (-1,1)
    %could just ignore reality constraints for Jac, and put A itself rand...
    for s=1:S(ell+1) %m=0
        A{ell+1}(s,ell+1)=(i^ell)*realAlpha{ell+1}(s,ell+1);
    end   
    for m=1:ell%m=/=0
        for s=1:S(ell+1)
          A{ell+1}(s,m+ell+1)=realAlpha{ell+1}(s,m+ell+1)+i*realAlpha{ell+1}(s,-m+ell+1);%m>0
          A{ell+1}(s,-m+ell+1)=((-1)^(ell+m))*(realAlpha{ell+1}(s,m+ell+1)-i*realAlpha{ell+1}(s,-m+ell+1));%m<0
        end
    end
end%got random volume coeffs in correct real subspace (tho reality not needed for Jac)

B = cell(P+1,1);
realBeta = cell(P+1,1);
for p=0:P
    B{p+1} = zeros(2*p+1,2*p+1);
    realBeta{p+1} = 2*rand(2*p+1,2*p+1)-ones(2*p+1,2*p+1);
    %could just ignore reality constraints for Jac, and put B itself rand...
    %I am ignoring positivity constraints (also OK for Jac, since positivity holds on a Euclidean open ball)
    B{p+1}(p+1,p+1)=realBeta{p+1}(p+1,p+1); %(u,v)=(0,0)
    for u = -p:p
        for v = -p:p
            if (u>0) || ((u==0) && (v>0))%(u,v)>=(0,0) in lex
                B{p+1}(u+p+1,v+p+1)=realBeta{p+1}(u+p+1,v+p+1)+i*realBeta{p+1}(-u+p+1,-v+p+1);
                B{p+1}(-u+p+1,-v+p+1)=((-1)^(u+v))*(realBeta{p+1}(u+p+1,v+p+1)-i*realBeta{p+1}(-u+p+1,-v+p+1));
            end
        end
    end  
end%got random distribution coeffs in correct real subspace (tho reality not needed for Jac)

%now try uniform in-plane
% for p=0:P
%     for u=-p:p
%         for v=-p:p
%             if v~=0
%                 B{p+1}(u+p+1,v+p+1)=0;
%             end
%         end
%     end
% end


curlyB = cell(2*L+1,2*L+1,L+1,L+1);%q1,q2,ell1,ell2
for q1=-L:L
    for q2=-L:L
        for ell1 = abs(q1):2:L %enforces ell1 >= |q1|, and ell1=q1(mod2)
            for ell2 = abs(q2):2:L %ditto
                if (max(abs(q1+q2),abs(ell1-ell2)) <= P)%enforces other bounds
                    %curlyB{q1+L+1,q2+L+1,ell1+1,ell2+1}=zeros(2*ell1+1,2*ell2+1);%indexed by m1,m2
                    M=zeros(2*ell1+1,2*ell2+1);%shorter notation, going to become curlyB cell above
                    for m1=-ell1:ell1
                        for m2=-ell2:ell2
                            for p=max([abs(m1+m2) abs(q1+q2) abs(ell1-ell2)]):min(ell1+ell2,P)
                               M(m1+ell1+1,m2+ell2+1)=M(m1+ell1+1,m2+ell2+1)+(B{p+1}(-m1-m2+p+1,-q1-q2+p+1))*(clebschgordan(ell1,ell2,p,m1,m2,m1+m2))*(clebschgordan(ell1,ell2,p,q1,q2,q1+q2))*((-1)^(m1+m2))*(1/(2*p+1));           
                            end
                        end
                    end
                    curlyB{q1+L+1,q2+L+1,ell1+1,ell2+1} = M;
                end
            end
        end 
    end 
end%got my curly B's for second moment....hopefully correct

curlylittleb=cell(2*L+1,L+1);%column vectors for first moment, indexed by q,ell
for q=-min(L,P):min(L,P)%changed from -L:L here!
    for ell=abs(q):min(L,P)
        curlylittleb{q+L+1,ell+1}=zeros(2*ell+1,1);
        for m=-ell:ell
            curlylittleb{q+L+1,ell+1}(m+ell+1,1)= (B{ell+1}(-m+ell+1,-q+ell+1))*((-1)^m)*(1/(2*ell+1));
        end
    end
end%got my curly b's for first moment....hopefully correct

%btw, no conjugates needed in my formulation here, following paper with no
%conjugates

%build second moment, volume block of Jacobian now
Jac_secondMoment_volA = zeros(sizeM2,sizeA);
columnCounter = 1;
%Jac_secondMoment_volA = [];
for ell=0:L
    for m=-ell:ell
        for s=1:S(ell+1)
            %col =[];%going to be ell,m,s column of this Jacobian block
            col = zeros(sizeM2,1);%going to be ell,m,s column of this Jacobian block
            rowCounter = 1;
            for q1=-L:L
                for q2=-L:q1%wlog q1>=q2, other moments exact repeats
                    if (abs(q1+q2)<=P)%otherwise moment identically 0
                        prod = zeros(T(abs(q1)+1),T(abs(q2)+1));
                        for ell1=abs(q1):2:L
                            for ell2=abs(q2):2:L
                                if (abs(ell1-ell2)<=P)
                                    if ell==ell1
                                        E_lms=zeros(S(ell+1),2*ell+1);
                                        E_lms(s,m+ell+1)=1;%1 in lms spot, 0 else
                                        prod=prod+Gamma{abs(q1)+1,ell1+1}*(E_lms)*curlyB{q1+L+1,q2+L+1,ell1+1,ell2+1}*transpose(A{ell2+1})*transpose(Gamma{abs(q2)+1,ell2+1});
                                    end
                                    if ell==ell2
                                        E_lms=zeros(S(ell+1),2*ell+1);
                                        E_lms(s,m+ell+1)=1;%1 in lms spot, 0 else
                                        prod=prod+Gamma{abs(q1)+1,ell1+1}*A{ell1+1}*curlyB{q1+L+1,q2+L+1,ell1+1,ell2+1}*transpose(E_lms)*transpose(Gamma{abs(q2)+1,ell2+1});
                                    end
                                    % GAMMAS HERE ARE WRONG Gamma{q1+L+1,ell1+1}*A{ell1+1}*curlyB{q1+L+1,q2+L+1,ell1+1,ell2+1}*transpose(A{ell2+1})*transpose(Gamma{q2+L+1,ell2+1})
                                    %keep at it!
                                end
                            end
                        end
                        if q1 > q2
                            prodVec = reshape(prod,T(abs(q1)+1)*T(abs(q2)+1),1);
                            col(rowCounter:rowCounter+size(prodVec,1)-1) = prodVec;
                            rowCounter = rowCounter + size(prodVec,1);
                            %col = [col; reshape(prod,T(abs(q1)+1)*T(abs(q2)+1),1)];
                        else %q1=q2, and should avoid duplicates by taking t1>=t2
                            mask = tril(true(T(abs(q1)+1)));
                            prodVec = prod(mask);
                            col(rowCounter:rowCounter+size(prodVec,1)-1) = prodVec;
                            rowCounter = rowCounter + size(prodVec,1);
                            %col = [col; prod(mask)];
                        end
                    end
                end
            end
            Jac_secondMoment_volA(:,columnCounter) = col;
            columnCounter = columnCounter + 1;
            %Jac_secondMoment_volA = [Jac_secondMoment_volA col];%got col now
        end
    end
end

%Jac_firstMoment_volA = [];
Jac_firstMoment_volA = zeros(sizeM1,sizeA);
columnCounter = 1;
for ell=0:L
    for m=-ell:ell
        for s=1:S(ell+1)
            %col =[];%going to be ell,m,s column of this Jacobian block
            col = zeros(sizeM1,1);%going to be ell,m,s column of this Jacobian block
            rowCounter = 1;
            for q1=-min(L,P):min(L,P)%changed from -L:L here!
                prod=zeros(T(abs(q1)+1),1);
                for ell1=abs(q1):2:min(L,P)%changed right limit from L here!
                    if ell==ell1
                        E_lms=zeros(S(ell+1),2*ell+1);
                        E_lms(s,m+ell+1)=1;%1 in lms spot, 0 else
                        prod=prod+Gamma{abs(q1)+1,ell1+1}*(E_lms)*curlylittleb{q1+L+1,ell1+1};
                    end
                end
                %col=[col;prod];
                col(rowCounter:rowCounter+size(prod,1)-1)=prod;
                rowCounter = rowCounter + size(prod,1);
            end
            %Jac_firstMoment_volA = [Jac_firstMoment_volA col];
            Jac_firstMoment_volA(:,columnCounter) = col;
            columnCounter = columnCounter + 1;
        end
    end
end

%which columns to keep for distB in uniform in-plane case
% keepers=[];
% counter=1;
% for p=0:P%just for uniform in-plane
%     for u=-p:p
%         for v=-p:p
%             if v==0
%                 keepers=[keepers counter];
%             end
%             counter=counter+1;
%         end
%     end
% end


%Jac_secondMoment_distB = [];
Jac_secondMoment_distB = zeros(sizeM2,sizeB);
columnCounter = 1;
for p=0:P
    for u=-p:p
        for v=-p:p
            %col=[];%going to be p,u,v column of this Jacobian block
            col = zeros(sizeM2,1);%going to be p,u,v column of this Jacobian block
            rowCounter = 1;
            for q1=-L:L
                for q2=-L:q1%wlog q1>=q2, other moments exact repeats, AS ABOVE (important to match row labels and order)
                    if (abs(q1+q2)<=P)%otherwise moment identically 0
                        prod = zeros(T(abs(q1)+1),T(abs(q2)+1));
                        if ((-q1-q2)==v)%otherwise no derivate wrt B_(p,u,v) here
                            for ell1=abs(q1):2:L
                                for ell2=abs(q2):2:L
                                    if (abs(ell1-ell2)<=P)
                                        middle=zeros(2*ell1+1,2*ell2+1);
                                        for m1=-ell1:ell1
                                            for m2=-ell2:ell2
                                                if ((-m1-m2)==u)%only these entries in middle have non0 deriv wrt B_(p,u,v)
                                                    if (max([abs(m1+m2) abs(q1+q2) abs(ell1-ell2)])<=p) && (p<=min(ell1+ell2,P))
                                                        middle(m1+ell1+1,m2+ell2+1)=(clebschgordan(ell1,ell2,p,m1,m2,m1+m2))*(clebschgordan(ell1,ell2,p,q1,q2,q1+q2))*((-1)^(m1+m2))*(1/(2*p+1));
                                                    end
                                                end
                                            end
                                        end
                                        %now have middle
                                        prod=prod+Gamma{abs(q1)+1,ell1+1}*A{ell1+1}*middle*transpose(A{ell2+1})*transpose(Gamma{abs(q2)+1,ell2+1});
                                    end
                                end
                            end
                        end
                        %now have prod (summed over l1,l2)
                        if q1 > q2
                            %col=[col;reshape(prod,T(abs(q1)+1)*T(abs(q2)+1),1)];
                            prodVec = reshape(prod,T(abs(q1)+1)*T(abs(q2)+1),1);
                            col(rowCounter:rowCounter+size(prodVec,1)-1) = prodVec;
                            rowCounter = rowCounter + size(prodVec,1);
                        else %q1=q2
                            mask = tril(true(T(abs(q1)+1)));
                            prodVec = prod(mask);
                            col(rowCounter:rowCounter+size(prodVec,1)-1) = prodVec;
                            rowCounter = rowCounter + size(prodVec,1);  
                        end
                    end
                end
            end
            %now have col
            %Jac_secondMoment_distB = [Jac_secondMoment_distB col];
            Jac_secondMoment_distB(:,columnCounter) = col;
            columnCounter = columnCounter + 1;
        end 
    end
end

%Jac_secondMoment_distB = Jac_secondMoment_distB([1:size(Jac_secondMoment_distB,1)],keepers);
%just for uniform in-plane

%Jac_firstMoment_distB = [];
Jac_firstMoment_distB = zeros(sizeM1,sizeB);
columnCounter = 1;
for p=0:P
    for u=-p:p
        for v=-p:p
            %col=[];%going to be p,u,v column of this Jacobian block
            col = zeros(sizeM1,1);%going to be p,u,v column of this Jacobian block
            rowCounter = 1;
            for q1=-min(L,P):min(L,P)
                prod=zeros(T(abs(q1)+1),1);
                if ((-q1)==v)%otherwise deriv wrt B_puv is 0
                    if ((p<=L)&&(p>=abs(q1))&&((mod(p-q1,2))==0))
                        rightVec = zeros(2*p+1,1);
                        rightVec(-u+p+1)=((-1)^u)*(1/(2*p+1));
                        prod=Gamma{abs(q1)+1,p+1}*A{p+1}*rightVec;
                    end
                end
                %now got prod
                %col=[col;prod];
                col(rowCounter:rowCounter+size(prod,1)-1)=prod;
                rowCounter = rowCounter + size(prod,1);
            end
            %now have col
            %Jac_firstMoment_distB = [Jac_firstMoment_distB col];
            Jac_firstMoment_distB(:,columnCounter)=col;
            columnCounter = columnCounter+1;
        end
    end
end

%Jac_firstMoment_distB = Jac_firstMoment_distB([1:size(Jac_firstMoment_distB,1)],keepers);
%just for uniform in-plane


Jac_total = [Jac_secondMoment_volA Jac_secondMoment_distB;
             Jac_firstMoment_volA Jac_firstMoment_distB];

singVal_Jac_total = svd(Jac_total)

myCond = singVal_Jac_total(1)/singVal_Jac_total(size(singVal_Jac_total,1)-3)

[~,~,V] = svd(Jac_total);
kerA = V(1:sizeA,end-2:end);%this is wrong if corank > 3
[U_A, ~, ~] = svd(kerA);
tanSpA = U_A(:,4:end);
tanSpAB = [tanSpA; zeros(sizeB,size(tanSpA,2))];
Proj = tanSpAB*transpose(tanSpAB);
invA = Proj*pinv(Jac_total);
singA = svd(invA);
myCondA = singA(1)/singA(sizeA-3)
%tanSp = V(1:sizeA,1:size(V,2)-3);




%runMin = 1000;
%for i=1:size(Jac_total,1)
%    for j=i+1:size(Jac_total,1)
%        sing = svd([Jac_total(i,:); Jac_total(j,:)]);
%        runMin = min(runMin, sing(2));
%    end
%end






















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a=A{1};
% maxL = size(a,2);
% for j=1:maxL
%     a{j} = rand(size(a{j})) + i*rand(size(a{j}));%the leftmost columns in a{j} should be purely real or purely imaginary depending on j
% end
% A{1}=a;
% maxP = size(B,1);
% 
% for j = 1:maxP
%     B{j} = rand(size(B{j})) + i*rand(size(B{j}));%messing up reality constraints on B, but shouldnt affect the Jacobian test
% end
% 
% tt0 = toc;
% fprintf('DONE in about %d seconds \n', round(tt0));
% 
% % main parameters
% L = max(gamma.band_idx_3d);  % the overall degree. Common to both 2D and 3D
% P = size(B,1);               % expansion length of the distribution
% M = max(gamma.ang_idx_2d)+1; % angular size of the moment
% vec_B = FromCellArr2Vec({1},B);
% vec_A = A_B_to_VecAB(A, [], M-1);
% 
% fprintf('*****  MAIN PARAMETERS: ****** \n');
% fprintf('Expansion length of the distribution: %d \n',P)
% fprintf('The overall degree (2D and 3D): %d \n',L);
% fprintf('\n \n');
% 
% lengthB = length(vec_B);
% lengthA = length(vec_A);
% 
% % preprocessing
% tic
% fprintf('Preprocessing..');
% [sign_mat, Gamma_mat] = preprocessing_mu1_coef(P, L, gamma);
% [C_tensor]            = preprocessing_mu2_coef_V2(gamma, P, M);
% [C_array] = mu2_C_coefs_PSWF_v1_fixed(L, P, M); % NO NEED. just for comparison with older version
% tt = toc;
% fprintf('DONE in about %d seconds \n', round(tt));
% 
% % calculate moments
% tic
% fprintf('Calculating moments..');
% mu1 = FirstMoment_PSWF_v2(A, B, Gamma_mat, sign_mat);
% mu2 = SecondMoment_PSWF_v2(A, B, gamma, C_tensor, Gamma_mat);
% tt2 = toc;
% fprintf('DONE in about %d seconds \n', round(tt2));
% 
% % now getting the coefs
% tic
% fprintf('Computing coefficients..');
% B_coef_mu1     = get_B_coefs_mu1(A, lengthB, P, Gamma_mat, sign_mat);
% B_coefs_tensor = getting_B_coefs_PSWF_V3(A, P, lengthB, gamma, C_array);
% [A_coef_mu1, conj_A_coef_mu1]           = get_A_coefs_mu1_full_A_vec(A, ...
%                                                 B, Gamma_mat, sign_mat);
% [A_coef_mu2, A_coef_mu2_conj, conj_ind] = get_A_coefs_mu2_full_vec(A,...
%                                                 B, Gamma_mat, C_tensor);
% tt3 = toc;
% fprintf('DONE in about %d seconds \n', round(tt3));
% 
% % getting into a matrix shape
% B_sys_mu1 = reshape(B_coef_mu1, numel(mu1),lengthB);
% B_sys_mu2 = reshape(B_coefs_tensor, numel(mu2),lengthB);
% 
% A_sys_mu1      = reshape(A_coef_mu1, numel(mu1),lengthA);
% A_sys_mu1_conj = reshape(conj_A_coef_mu1, numel(mu1),lengthA);
% A_sys_mu2      = reshape(A_coef_mu2, numel(mu2),lengthA);
% A_sys_mu2_conj = reshape(A_coef_mu2_conj, numel(mu2),lengthA);
% 
% A_sys_mu1_conj = A_sys_mu1_conj(:,conj_ind);
% A_sys_mu2_conj = A_sys_mu2_conj(:,conj_ind);
% 
% 
% % debugging the linear parts, 
% test_lin_parts = 0;
% if test_lin_parts
%     norm(B_sys_mu1*vec_B - mu1(:))
%     norm(B_sys_mu2*vec_B - mu2(:))
%     norm(A_sys_mu1*vec_A + A_sys_mu1_conj*conj(vec_A) - mu1(:))
% end
% 
% %% getting the parameter of each ecolumn in the system
% 
% % A system
% A_l_deg = A;
% A_s_deg = A;
% A_n_deg = A;
% for j=1:length(A{1})
%     A_l_deg{1}{j} = ones(size(A_l_deg{1}{j}))*(j-1);
%     A_s_deg{1}{j} = diag(1:size(A_l_deg{1}{j},1))*ones(size(A_l_deg{1}{j}));
%     A_n_deg{1}{j} = ones(size(A_l_deg{1}{j}))*diag(1:size(A_l_deg{1}{j},2));
% end
% l_deg_vecA = A_B_to_VecAB(A_l_deg, [], M-1);
% s_deg_vecA = A_B_to_VecAB(A_s_deg, [], M-1);
% n_deg_vecA = A_B_to_VecAB(A_n_deg, [], M-1);
% 
% % including the negative parameters
% l_deg_vecA = [l_deg_vecA(conj_ind); l_deg_vecA];
% s_deg_vecA = [s_deg_vecA(conj_ind); s_deg_vecA];
% n_deg_vecA = [-n_deg_vecA(conj_ind); n_deg_vecA];
% 
% % B system
% B_deg = B;
% B_u = B;
% B_v = B;
% for j=1:P
%     B_deg{j} = ones(size(B{j}))*(j-1);
%     B_u{j} = diag(1:size(B{j},1))*ones(size(B{j}));
%     B_v{j} = ones(size(B{j}))*diag(1:size(B{j},2));
% end
% deg_vecB = FromCellArr2Vec({0},B_deg);
% u_vecB = FromCellArr2Vec({1},B_u);
% v_vecB = FromCellArr2Vec({1},B_v);
% 
% 
% %% testing the systems and report
% B_system = [B_sys_mu2 ;B_sys_mu1];
% 
% A_system = [A_sys_mu1_conj, A_sys_mu1; ...
%             A_sys_mu2_conj , A_sys_mu2];
% 
% Jacobian_mat = [A_system, B_system];
% 
% r_b = rank(B_system);
% c_b = cond(B_system);
% 
% r_a = rank(A_system);
% c_a = cond(A_system);
% 
% [u,s,v] = svd([A_system]);
% singular_vals = diag(s);
% %hist(singular_vals(end-5:end),10); title('Least singular values -- volume');
% figure; scatter(1:numel(singular_vals), log10(singular_vals),'filled');
% title('Singular values -- volume');
% 
% [u2,s2,v2] = svd(Jacobian_mat);
% singular_vals2 = diag(s2);
% %hist(singular_vals2(end-3:end),10); title('Least singular values -- both vol and dist.');
% 
% 
% figure; scatter(1:numel(singular_vals2), log10(singular_vals2),'filled');
% title('Singular values -- BOTH vol. and dist.');
% 
% fprintf('----------------------------------------------------- \n');
% fprintf('We have %d degrees of freedom due to conjugation, and %d zero singular values \n', ...
%     numel(conj_ind), nnz(singular_vals2<1e-15));
% fprintf('----------------------------------------------------- \n');
% 
% 
% 
% fprintf('\n \n');
% fprintf('----------------------------------------------------- \n');
% fprintf('RESULTS:\n');
% fprintf('The distribution part, %d unknowns, (rank, condition): (%d, %e) \n' , lengthB, r_b, c_b)
% fprintf('The volume part, %d unknowns, (rank, condition): %d, %e \n',lengthA+numel(conj_ind), r_a, c_a)
% fprintf('----------------------------------------------------- \n');
% 


