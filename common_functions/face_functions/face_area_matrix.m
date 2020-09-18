%% Computes matrix M with faces' area on diagonal

%INPUTS:
    %V = vertices
    %F = faces
%OUTPUTS:   
    %M = 3|F|-by-3|F| sparse area matrix , M(i,i) = 1/3 * area(face_i)
    
%example usage:  ||u||^2 = u^t * M * u (mass-norm)
    %note: M is a positive definite diagonal matrix for all non-degenerate
    %meshes
    
function M = face_area_matrix(V,F)   
    dblA = doublearea(V,F);
    N = 3*length(F);
    M = sparse(1:N, 1:N, 1/6*repelem(dblA,3));
end