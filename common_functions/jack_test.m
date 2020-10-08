%jack's testing script, it might get moved to my folder later
%works with normal old face functions
mkdir('vertex_functions')
addpath('vertex_functions')
mkdir('face_functions')
addpath('face_functions')
x = 0:.1:1;
y = x';
[x,y] = meshgrid(x,y);
x = x(:);
y = y(:);
F = delaunay(x,y);
V = [x y];
max_iterations = 5;
[m,~] = size(V);
%u_prev = zeros(m*(max_iterations - 1), 1);
%for i = 1:max_iterations
%    u = laplace_eq_2D_quadprog_iterations(V, F, u_prev, i);
%    u_prev(i*(m-1) + 1: m*i) = u;
%end
%u1 = laplace_eq_2D_seg_quadprog(V, F);
%face_plotting(V,F,u1)
u = laplace_disc(V,F);
face_plotting(V,F,u)


function u = laplace_eq_2D_quadprog_iterations(V, F, u_prev, iterations)
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
    
    [n,~] = size(V);
    [k, ~] = size(F);
    err = 0.1;
    if iterations == 1
        G = [grad(V,F) speye(2*k, (iterations-1)*n)];
        M = speye(2*k);
        Mu = speye(n,n);
        c = rand(n,1);
        c = c/norm(c);
        u = quadprog(2*G'*M*G, [], [], [], trans(M*c), ones(n,1));
        while abs(u'*M*u - ones(n,1)) >= err
            c = u/norm(u);
            u = quadprog(2*G'*M*G, [], [], [], trans(M*c), ones(n,1));
        end
    end
    
end

function u = laplace_disc(V, F)
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
   
    L = -face_build_discontinuity_laplacian(V,F); %|F| * 3 x |F| * 3
    O_right = zeros(3*k,2*m); 
    O_bottom = zeros(2*m,3*k+2*m);
    G_tilde = [L O_right; O_bottom ];
    
    % generate discontinuity matrix
    V = [V zeros(length(V))];
    D = discontinuity(V,F); %D = 2|E|-by-3|F| sparse matrix , E = #edges
    V = V(:,1:2);
    
    % generate vector f used in the constraint
    u_placeholder = zeros(3*k,1);
    t_placeholder = ones(2*m,1);
    f = [ u_placeholder; t_placeholder ];
    
    % generate |E|x|E| identity matrix used in the constraint block matrix
    I = speye(2*m);
    
    % generate inequality constraint matrix and vector
    A = [D -I; -D -I];
    b = zeros(4*m,1);
    
    %boundary conditions
    B = unique(reshape(outline(F),[],1));
    BP1 = B(V(B,1) == 0);
    BP2 = B(V(B,1) == 1);
    
    %bad way to get locations of boundary vertices in face based function
    %synatx fix later
    BP1new = [];
    for i = 1:length(BP1)
        [r,s] = size(F);
        for j = 1:r
            for k = 1:s
                if F(j,k) == BP1(i)
                    BP1new = [BP1new 3*(j-1) + k];
                end
            end
        end
        
    end
    BP1new = BP1new';
    BP2new = [];
    for i = 1:length(BP2)
        [r,s] = size(F);
        for j = 1:r
            for k = 1:s
                if F(j,k) == BP2(i)
                    BP2new = [BP2new 3*(j-1) + k];
                end
            end
        end
        
    end
    BP2new = BP2new';
    [k,~] = size(F);
    
    equalityConstraint = sparse(BP2new,ones(length(BP2new),1),ones(length(BP2new),1), 3*k + 2*m, 1);
    equalityMatrix = sparse([BP1new; BP2new], [BP1new; BP2new], ones(length(BP1new) + length(BP2new), 1), 3*k + 2*m, 3*k + 2*m);
    
    y = quadprog(G_tilde,0.5*f,A,b,equalityMatrix,equalityConstraint);
    
    u = y(1:3*k);
    %Constructing a function for debugging
    v = ones(length(u),1);
    for i = 1:length(BP1)
        [r,~] = size(F);
        for j = 1:r
            if (BP1(i) == F(j,1) || BP1(i) == F(j,2) || BP1(i) == F(j,3))
                v(3*(j-1) + 1) = 0;
                v(3*(j-1) + 2) = 0;
                v(3*(j-1) + 3) = 0;
            end
        end
    end
norm(L*v,2) + norm(D*v,1)
norm(L*u,2) + norm(D*u,1)
end