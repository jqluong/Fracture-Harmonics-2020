function u = laplace_eq_2D(V, F, B)
    % solve laplace's equation on a 2D tiangle mesh surface
    %
    % V is vertex list
    %
    % F is face list
    %
    % B is boundary conditions. B should be a k by 3 matrix, where the
    % first column is x coordinate, second column is y coordinate, and
    % third column is the value of u at the corresponding (x, y)
    %
    % output is the solution u of \Delta u = 0 given boundary condition B
    
    % extract boundary condition
    g = B(:,3);
    
    % find indices for boundary condition vertices
    [~, V_bpos] = ismember(B(:,1:2), V, 'rows');
    
    % generate cotangent matrix
    % notation: L = cotangent matrix, subscript b for boundary (corresponds
    % to \partial in notes), and subscript i for interior (corresponds to
    % \Omega in notes)
    L = cotmatrix(V,F);
    % find column indices corresponding to non-boundary condition vertices
    L_ip = setdiff(1:length(L),V_bpos);
    % find submatrix that operates only on non-boundary condition vertices
    L_ii = L(L_ip, L_ip);
    
    u = zeros(length(V),1);
    u(V_bpos) = g;

    % if u_int = 0, L*u will have contributions from L_ib*u_b and L_bb*u_b,
    % so perform multiplication and extract values corresponding to 
    % non-boundary condition points (L_ib*u_b)
    b = L*u;
    b = b(L_ip);
    
    u_int = -1 * L_ii \ b;
    
    u(setdiff(1:end,V_bpos)) = u_int;
end