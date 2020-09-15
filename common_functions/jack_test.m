%jack's testing script, it might get moved to my folder later
%works with normal old face functions
mkdir('vertex_functions')
addpath('vertex_functions')
x = 0:.1:1;
y = x';
[x,y] = meshgrid(x,y);
x = x(:);
y = y(:);
F = delaunay(x,y);
V = [x y];
max_iterations = 5;
[m,~] = size(V);
u_prev = zeros(m*(max_iterations - 1), 1);
for i = 1:max_iterations
    u = laplace_eq_2D_quadprog_iterations(V, F, u_prev, i);
    u_prev(i*(m-1) + 1: m*i) = u;
end
%u = laplace_eq_2D_seg_quadprog(V,F);
%face_plotting(V,F,u)


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