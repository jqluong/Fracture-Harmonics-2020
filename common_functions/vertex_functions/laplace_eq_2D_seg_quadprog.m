function u = laplace_eq_2D_seg_quadprog(V, F)
    % solve laplace's equation on a 2D triangle mesh surface using quadprog
    % for 'segmented' functions, i.e. functions with face edges that may
    % not match
    %
    % V is matrix of vertex position
    %
    % F is matrix of face indices
    %
    % B is boundary conditions. B should be a k by 2 matrix, where the
    % first column is vertex position in V, and
    % second column is the value of u at the corresponding (x, y)
    %
    % output is the solution u of \Delta u = 0 given boundary condition B
    %
    % the format of solution u is u = [ u(F1_v1) u(F1_v2) u(F1_v3) u(F2_v1) u(F2_v2) u(F2_v3) ... ]
    
    E = edges(F); % edges matrix
    
    [m,~] = size(E); 
    
    [k,~] = size(F);
   
    % generate gradient matrix to act on [ u t ]
    G = face_grad(V,F);  % size(G) = 2k x 3k
    O = zeros(2*k);
    G_tilde = [ G O ];
    
    % generate discontinuity matrix
    D = face_discontinuity_matrix(V,F);
    
    % generate vector f used in the constraint
    u_placeholder = zeros(3*k,1);
    t_placeholder = ones(m,1);
    f = [ u_placeholder t_placeholder ];
    
    % generate |E|x|E| identity matrix used in the constraint block matrix
    I = speye(m);
    
    % generate inequality constraint matrix and vector
    A = [ -D -I; D I];
    b = zeros(2*m,1);
    
    y = quadprog(G_tilde,f,A,b);
    
    u = y(1:3*m);
    
    
end