function E = face_dirichlet_energy(V, F, u)
    % compute the Dirichlet Energy of a face defined function u
    % inputs: Mesh V, F, and function u
    % u is represented as:
    %   u = [ u(F1_v1) u(F1_v2) u(F1_v3) u(F2_v1) u(F2_v2) u(F2_v3) ... ]'
    % output is the energy
    
    A = doublearea(V, F);
    du = face_grad(V, F)*u;
    du = reshape(du, [], size(V,2)); % reshape du from dim*F x 1 to F x dim
    du = sum(du.^2, 2); % norm square of each gradient vector
    
    E = sum(A.*du)/4;
end