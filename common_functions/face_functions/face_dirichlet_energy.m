function E = face_dirichlet_energy(V, F, u, G)
    % compute the Dirichlet Energy of a face defined function u
    % inputs: Mesh V, F, function u, optional gradient G
    % u is represented as:
    %   u = [ u(F1_v1) u(F1_v2) u(F1_v3) u(F2_v1) u(F2_v2) u(F2_v3) ... ]'
    % output is the energy
    
    if nargin == 3 % gradient matrix not provided
        du = face_grad(V, F)*u;
    elseif nargin == 4 % gradient matrix provided
        du = G*u;
    else
        error('check all input arguments provided');
    end
    
    A = doublearea(V, F);
    du = reshape(du, size(V,2), [])'; % reshape du from dim*F x 1 to F x dim
    du = sum(du.^2, 2); % norm square of each gradient vector
    
    E = sum(A.*du)/4;
end