%TOTALLY NON-UNIFORM NOW

L=7; % put herewhatever is the volume degree in the PSWF expansion
P=L; %we're going to square these linear combs of Wigner entries up to this P
%so really this is "half P"

NumberOfSquares = 2;


B = cell(2*P+1,1);
for p = 0:2*P
    B{p+1} = zeros(2*p+1,2*p+1);
end


for t = 1:NumberOfSquares
    linComb = cell(P+1,1); %thing we're going to square (or rather, multiply by its conjugate)
    for p = 0:P
        linComb{p+1} = (2*rand(2*p+1,2*p+1)-ones(2*p+1,2*p+1)) + 1i*(2*rand(2*p+1,2*p+1)-ones(2*p+1,2*p+1));
    end
    for p1=0:P
        for p2=0:P
            for u1=-p1:p1
                for v1=-p1:p1
                    for u2=-p2:p2
                        for v2=-p2:p2
                            for p3=max([abs(p1-p2) abs(u1-u2) abs(v1-v2)]):(p1+p2)
                                B{p3+1}(u1-u2+p3+1,v1-v2+p3+1)=B{p3+1}(u1-u2+p3+1,v1-v2+p3+1)+linComb{p1+1}(u1+p1+1,v1+p1+1)*conj(linComb{p2+1}(u2+p2+1,v2+p2+1))*((-1)^(u2+v2))*clebschgordan(p1,p2,p3,u1,-u2,u1-u2)*clebschgordan(p1,p2,p3,v1,-v2,v1-v2); 
                            end
                        end
                    end
                end
            end
        end
    end
end

                            
            
                    
normalization = B{1}(1);
%B=cell(2*P+1,1);
for p=0:2*P
    %B{p+1} = zeros(2*p+1,2*p+1);
    B{p+1} = (1/normalization)*B{p+1};
end

%Bsize = zeros(2*P+1,1);
%for p=0:2*P
%    Bsize(p+1) = norm(B{p+1},'fro');
%end

save('B15.mat', 'B');
B
    