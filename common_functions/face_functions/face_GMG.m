%% Computes matrix GMG = G^t * M * G (for segmented functions)
    %u^t * GMG * u  = (||grad_u||_{L2})^2 = 
    %Dirichlet energy of u = 1/2 * u^t * GMG * u
    
%INPUTS
    %V = vertices
    %F = faces
%OUTPUTS
    %GMG = 3|F|-by-3|F| sparse matrix
    
    function [GMG] = face_GMG(V,F)
    
    G = face_grad(V,F);
    M = 3*face_area_matrix(V,F);
    GMG = sparse((G')*M*G);
    
    end