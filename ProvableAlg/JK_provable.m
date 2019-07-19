% IMPORTANT PSWF PARAMETERS MATCH WHATEVEVER IS USED TO GENERATE A !

% parameters
gridSize = 23;   % volume resolution
delta    = .99;
beta     = 1; % between 0 and 1,


gridSize = 43;   % volume resolution
delta    = .99;
beta     = 0.82; % between 0 and 1, controls the expansion length

% getting the map -- volume expansion coefficients
map    = ReadMRC('emd_0409.mrc'); % upload the 128X128X128 vol
vol    = cryo_downsample(map,[gridSize gridSize gridSize]); % Downsamples map
eps_p  = 1e-3;    % Prescribed accuracy in images side
radius = floor(gridSize/2);
c      = beta*pi*radius;              % Nyquist bandlimit if beta == 1
A = pswf_t_f_3d(vol, beta, delta); %previously call A_full but now truncation is below
gamma  = PSWF_2D_3D_T_mat(c, delta, eps_p);

numA = 0;
for ell=0:L
    numA = numA + (2*ell+1)*S(ell+1);
end




% gridSize = 43; 
% radius = floor(gridSize/2);
% beta=.8;
% c = beta*pi*radius;
% delta    = .99;     
% eps_p = 1e-3;
% gamma = PSWF_2D_3D_T_mat(c, delta, eps_p);
% 
% load('volume_coefs_and_gamma','A') % MATCH WITH PSWF PARAMETERS ABOVE ! 
% % HERE A_ell IS IN NIR'S EFFICIENT FORMAT WITH ell+1 COLUMNS !


L=-1;
S=[];
while ((2*(L+1)+1) <= sum(gamma.band_idx_3d(:)==L+1))
    S = [S sum(gamma.band_idx_3d(:)==L+1)];
    L=L+1;
end
% S = zeros(1,L+1);
% for i=1:L+1 %really 0 to L
%     S(i)=sum(ind3D(:)==i-1);
% end
ind2D = gamma.ang_idx_2d;
T = zeros(1,L+1);
for i=1:L+1 %really 0 to L
    T(i)=sum(ind2D(:)==i-1);
end
Gamma=cell(L+1,L+1);%gamma^q,ell with q>=0.  Use gamma^q,ell = gamma^-q,ell..


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


nirA = A{1}';
ourA = cell(L+1,1);
for ell=0:L 
    nirAell = nirA{ell+1,1};
    ourAell = zeros(size(nirAell,1),2*ell+1);
    for col=0:ell
        ourAell(:,col+ell+1)=nirAell(:,col+1);
    end
    for col=-ell:-1
        ourAell(:,col+ell+1)=((-1)^(ell+col))*conj(nirAell(:,-col+1));
    end
    ourA{ell+1,1}=ourAell;
end
A = ourA;


load('B15') % MUST BE TOTALLY NON-UNIFORM !  MUST HAVE degree >= 2*L !!
P = size(B,1)-1;




% sizeA = 0;
% for ell = 1:L+1 %really 0:L 
%     sizeA = sizeA + (2*ell-1)*S(ell);
% end
% 
% indA = zeros(sizeA,3);
% counter=1;
% for ell=0:L 
%     for m=-ell:ell
%         for s=1:S(ell+1)
%             indA(counter,:)=[ell m s];
%             counter=counter+1;
%         end
%     end
% end
% 
% sizeB = 0;
% for p = 0:P
%     sizeB = sizeB + (2*p+1);
% end %count for totally non-uniform %BUT THIS IS KNOWN ANYHOW HERE
% 
% maxT = max(T);
% 
% sizeM1 = T(1);
% 
% sizeM2 = 0;
% for q1 = 0:L %q2 = -q1 because in-plane uniform, wlog q1 >= q2
%     if  q1 == 0
%         sizeM2 = sizeM2 + nchoosek(T(abs(q1)+1)+1,2);
%     else
%         sizeM2 = sizeM2 + T(abs(q1)+1)*T(abs(q1)+1);
%     end
% end
%     
    

% A = cell(L+1,1);
% for ell=0:L
%     A{ell+1}=zeros(S(ell+1),2*ell+1);
% end
% for ell=0:L
%     for m=0:ell
%         A{ell+1}(:,m+ell+1)=(2*rand(S(ell+1),1)-ones(S(ell+1),1)) + 1i*(2*rand(S(ell+1),1)-ones(S(ell+1),1));
%     end
%     for m=-ell:-1
%         A{ell+1}(:,m+ell+1) = ((-1)^(ell+m))*conj(A{ell+1}(:,-m+ell+1));
%     end
% end
    


% B = cell(P+1,1);
% realBeta = cell(P+1,1);
% for p=0:P
%     B{p+1} = zeros(2*p+1,2*p+1);
%     realBeta{p+1} = 2*rand(2*p+1,2*p+1)-ones(2*p+1,2*p+1);
%     %could just ignore reality constraints for Jac, and put B itself rand...
%     %I am ignoring positivity constraints (also OK for Jac, since positivity holds on a Euclidean open ball)
%     B{p+1}(p+1,p+1)=realBeta{p+1}(p+1,p+1); %(u,v)=(0,0)
%     for u = -p:p
%         for v = -p:p
%             if (u>0) || ((u==0) && (v>0))%(u,v)>=(0,0) in lex
%                 B{p+1}(u+p+1,v+p+1)=realBeta{p+1}(u+p+1,v+p+1)+1i*realBeta{p+1}(-u+p+1,-v+p+1);
%                 B{p+1}(-u+p+1,-v+p+1)=((-1)^(u+v))*(realBeta{p+1}(u+p+1,v+p+1)-1i*realBeta{p+1}(-u+p+1,-v+p+1));
%             end
%         end
%     end  
% end%got random distribution coeffs in correct real subspace (tho reality not needed for Jac)

%now try uniform in-plane %JUST ZEROS OUT IN CASE NOT IN-PLANE INVARIANT
%  for p=0:P
%      for u=-p:p
%          for v=-p:p
%              if v~=0
%                  B{p+1}(u+p+1,v+p+1)=0;
%              end
%          end
%      end
%  end


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



M2 = cell(2*L+1,2*L+1);
for q1 = -L:L
    for q2 = -L:L %here I don't assume q1>=q2 or <=
        M2{q1+L+1,q2+L+1} = zeros(T(abs(q1)+1),T(abs(q2)+1));
        for ell1 = abs(q1):2:L
            for ell2 = abs(q2):2:L
                if abs(ell1 - ell2) <= P
                    M2{q1+L+1,q2+L+1} = M2{q1+L+1,q2+L+1} + Gamma{abs(q1)+1,ell1+1}*A{ell1+1}*curlyB{q1+L+1,q2+L+1,ell1+1,ell2+1}*transpose(A{ell2+1})*transpose(Gamma{abs(q2)+1,ell2+1});
                end
            end
        end    
    end
end


% noise=0;
% for q1 = -L:L
%     for q2 = -L:L
%         M2{q1+L+1,q2+L+1} = M2{q1+L+1,q2+L+1} + noise*norm(M2{q1+L+1,q2+L+1},'fro')*rand(T(abs(q1)+1),T(abs(q2)+1));
%     end
% end


curlylittleb=cell(2*L+1,L+1);%column vectors for first moment, indexed by q,ell
for q=-min(L,P):min(L,P)%changed from -L:L here!
    for ell=abs(q):min(L,P)
        curlylittleb{q+L+1,ell+1}=zeros(2*ell+1,1);
        for m=-ell:ell
            curlylittleb{q+L+1,ell+1}(m+ell+1,1)= (B{ell+1}(-m+ell+1,-q+ell+1))*((-1)^m)*(1/(2*ell+1));
        end
    end
end%got my curly b's for first moment....hopefully correct


M1 = cell(2*L+1,1);
for q = -min(L,P):min(L,P) %note need to spell this out in paper still!
    M1{q+L+1}=zeros(T(abs(q)+1),1);
    for ell=abs(q):2:min(L,P)
        M1{q+L+1} = M1{q+L+1} + Gamma{abs(q)+1,ell+1}*A{ell+1}*curlylittleb{q+L+1,ell+1};
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic
myA = cell(L+1,1);
for ell=0:L
    myA{ell+1}=zeros(S(ell+1),2*ell+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B1 = curlyB{2*L+1,2*L+1,L+1,L+1};
B2 = curlyB{2*L+1,1,L+1,L+1};
[Q12,D12] = eig(B1*inv(B2));
b12 = (inv(Q12))*(curlylittleb{2*L+1,L+1});
[X12,E12] = eig((M2{2*L+1,2*L+1})*pinv(M2{2*L+1,1}));

ratios = diag(E12) * transpose(1./diag(D12));

permutation12 = (abs(ratios-1) < 10^(-5));

X12=X12*permutation12;

Lambda12 = diag((pinv(X12)*(M1{2*L+1}))./b12);

myA{L+1} = pinv(Gamma{L+1,L+1})*X12*Lambda12*inv(Q12);

for q1 = L-2:-2:0 %q2=L in recursive solve
    mySum = zeros(T(q1+1),T(L+1));
    for ell1 = q1+2:2:L %these are the "high ell1"
        mySum = mySum + Gamma{q1+1,ell1+1}*myA{ell1+1}*curlyB{q1+L+1,2*L+1,ell1+1,L+1}*transpose(myA{L+1})*transpose(Gamma{L+1,L+1});
    end
    myA{q1+1}=pinv(Gamma{q1+1,q1+1})*(M2{q1+L+1,2*L+1}-mySum)*pinv(curlyB{q1+L+1,2*L+1,q1+1,L+1}*transpose(myA{L+1})*transpose(Gamma{L+1,L+1}));
end
%got A_ell for all ell = L mod (2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B3 = curlyB{2*L,2*L,L,L};
B4 = curlyB{2*L,2,L,L};
[Q34,D34] = eig(B3*inv(B4));
b34 = (inv(Q34))*(curlylittleb{2*L,L});


[X34,E34] = eig((M2{2*L,2*L})*pinv(M2{2*L,2}));

ratios = diag(E34) * transpose(1./diag(D34));

permutation34 = (abs(ratios-1) < 10^(-5));

X34=X34*permutation34;

Lambda34 = diag((pinv(X34)*(M1{2*L}))./b34);

myA{L} = pinv(Gamma{L,L})*X34*Lambda34*inv(Q34);

for ell = L-3:-2:0 %change this for other parity !!! %should probably switch to better notation like above but w/e
    mySum = zeros(T(ell+1),T(L));
    for otherEll = ell+2:2:L-1
        mySum = mySum + Gamma{ell+1,otherEll+1}*myA{otherEll+1}*curlyB{ell+L+1,2*L,otherEll+1,L}*transpose(myA{L})*transpose(Gamma{L,L});
    end
    myA{ell+1}=pinv(Gamma{ell+1,ell+1})*(M2{ell+L+1,2*L}-mySum)*pinv(curlyB{ell+L+1,2*L,ell+1,L}*transpose(myA{L})*transpose(Gamma{L,L}));
end
toc
%now compare myA (computed coefficients) with A (ground truth)
relErrors = zeros(1,L+1);
for ell=0:L
    relErrors(ell+1) = norm(A{ell+1} - myA{ell+1}, 'fro')/norm(A{ell+1}, 'fro');
end
relErrors


totalAnorm = 0;
for ell=0:L
    totalAnorm = totalAnorm + norm(A{ell+1},'fro')^2;
end
totalAnorm = sqrt(totalAnorm);

totalDiffnorm = 0;
for ell=0:L
    totalDiffnorm = totalDiffnorm + norm(A{ell+1}-myA{ell+1},'fro')^2;
end
totalDiffnorm = sqrt(totalDiffnorm);

totRelErr = totalDiffnorm / totalAnorm


 %totalRelErr 5.4466e-11;
 %Elapsed time is 0.237955 seconds.
 %L=7, P =14
 %totalNumUnknowns = 1080
 %NumEqns??
 %
 
 totalNumUnknowns = 0;
 for ell=0:L
     totalNumUnknowns = totalNumUnknowns + (2*ell+1)*S(ell+1);
 end