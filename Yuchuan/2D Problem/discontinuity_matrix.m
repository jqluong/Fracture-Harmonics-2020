%% Computes discontinuity operator D 

%requires gptoolbox

%INPUTS:
    %F = faces, |F|-by-3 matrix 
%OUTPUTS:
    %D = 2|E|-by-3|F| sparse matrix , E = #edges
    
%D*u(i) and D*u(i+|E|) measures the discontinuity at edge i (at the two endpoints)
                                


function [D] = discontinuity_matrix(F)

temp = F';
f_vert = temp(:);
E = edges(F); %Edges
ET = edge_triangle_adjacency(F,E);
c = zeros(size(E,2),2);
d = zeros(size(E,2),2);

for i = 1:size(E,1)
    if ET(i,2) == -1
        c(i,:) = -1;
        d(i,:) = -1;
    else 
        for j = 0:2
            if f_vert(3*ET(i,1)-j)== E(i,1)
                c(i,1) = 3*ET(i,1)-j;
            elseif f_vert(3*ET(i,1)-j)== E(i,2)
                c(i,2) = 3*ET(i,1)-j;
            end
        end
        
        for k = 0:2
            if f_vert(3*ET(i,2)-k)== E(i,1)
                d(i,1) = 3*ET(i,2)-k;
            elseif f_vert(3*ET(i,2)-k)== E(i,2)
                d(i,2) = 3*ET(i,2)-k;
            end
        end 
    end
end
    
N = size(E,1);
D = zeros(2*N,length(f_vert));
    
for m = 1:N
    if c(m,1) <= 0
        continue;
    else 
        D(m,c(m,1)) = 1;
        D(m,d(m,1)) = -1;
        D(m+N,c(m,2)) = 1;
        D(m+N,d(m,2)) = -1;
    end
    
end
    
D = sparse(D);
    
end













