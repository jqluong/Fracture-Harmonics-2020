%% Computes weighted discontinuity operator D 

%requires gptoolbox

%INPUTS:
    %F = faces, |F|-by-3 matrix 
%OUTPUTS:
    %D = 2|E|-by-3|F| sparse matrix , E = #edges
    
%D*u(i) and D*u(i+|E|) measures the discontinuity at edge i (at the two endpoints)
                                
function [D] = discontinuity_matrix(V,F)
    temp = F';
    f_vert = temp(:);
    E = edges(F); %Edges
    ET = edge_triangle_adjacency(F,E);
    c = zeros(size(E,2),2);
    d = zeros(size(E,2),2);
    lE = length(E);
    lFv = length(F_vert);

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

    Dpos = zeros(2*lE, 3);
    ii = 1;
    for m = 1:lE
        if c(m,1) <= 0
            continue;
        else
            Dpos_c = [m c(m,1) 1; m d(m,1) -1; (m + lE) c(m,2) 1; (m + lE) d(m,2) -1];
            Dpos(ii:ii+3,:) = Dpos_c;
        end
        ii = ii+4;
    end
    ii = ii-1;
   
    W = face_edge_lengths_matrix(V,F);
    W = kron(speye(2), W);
    D = W*sparse(Dpos(1:ii,1), Dpos(1:ii,2), Dpos(1:ii,3), 2*lE, lFv);
end