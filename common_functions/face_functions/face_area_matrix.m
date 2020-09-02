%% Computes matrix M with faces' area on diagonal

%INPUTS:
    %V = vertices
    %F = faces
%OUTPUTS:   
    %M = 3|F|-by-3|F| sparse area matrix , M(i,i) = 1/3 * area(face_i)
    
%example usage:  ||u||^2 = u^t * M * u (mass-norm)
    %note: M is a positive definite diagonal matrix for all non-degenerate
    %meshes
    
    function [M] = face_area_matrix(V,F)
    
    dblA = (1/2)*doublearea(V,F);
    rep_area = repelem(dblA,3);
    M = sparse((1/3)*diag(rep_area));

    end