%% Computes a mass matrix for face functions

%INPUTS:
    %V = vertices
    %F = faces
%OUTPUTS:   
    %M = 3|F|-by-3|F| sparse mass matrix

function M = face_massmatrix(V,F)  
    A = (1/2)*doublearea(V,F);
    M = sparse(3*size(F,1),3*size(F,1));
    for i = 1:size(F,1)
        B = zeros(3,3);
        B(:) = (1/12)*A(i);
        M = M + kron(sparse(i,i,1,size(F,1),size(F,1)) , B);
    end
    
    C = sparse(1:3*size(F,1), 1:3*size(F,1), (1/12)*repelem(A,3));
    M = M + C;
end