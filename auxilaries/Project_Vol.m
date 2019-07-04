function [A_proj] = Project_Vol(A)
%

L = length(A{1});
A_proj = cell(1);
A_proj{1} = cell(1,L);
for j=1:L
    A_proj{1}{j} = A{1}{j};
    if mod(j,2)==0
        for l = 1:2:j
            A_proj{1}{j}(:,l) = 0;
        end
    else
        A_proj{1}{1} = real(A_proj{1}{1});
        for l = 2:2:j
            A_proj{1}{j}(:,l) = 0;
        end
    end
end
end

