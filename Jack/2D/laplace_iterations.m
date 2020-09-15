%jack's testing script, it might get moved to my folder later
%works with normal old face functions
mkdir('vertex_functions')
addpath('vertex_functions')
%x = 0:.1:1;
%y = x';
%[x,y] = meshgrid(x,y);
%x = x(:);
%y = y(:);
s = 100;
theta = [0:2*pi()/s:2*pi()]';
x = cos(theta);
y = sin(theta);
%[x,y] = meshgrid(x,y);
F = delaunay(x,y);
V = [x y];
max_iterations = 5;
[m,~] = size(V);
u_prev = zeros(m, max_iterations - 1);
for i = 1:max_iterations
    u = laplace_eq_2D_quadprog_iterations(V, F, u_prev, i, max_iterations);
    u_prev(:,i) = u;
end
%u = laplace_eq_2D_seg_quadprog(V,F);
%face_plotting(V,F,u)


function u = laplace_eq_2D_quadprog_iterations(V, F, u_prev, iterations, max_iterations)
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
    err = 10^-5;
    if iterations == 1
        G = grad(V,F);
        %M = speye(2*k);
        cellArea = repelem(power(doublearea(V, F)/2, 0.5), 2);
        M = spdiags(cellArea, 0, k*2, k*2);
        Mu = speye(n,n);
        c = rand(n,1);
        c = c/norm(c);
        u = quadprog(2*G'*M*G, [], [], [], transpose(Mu*c), 1);
        while abs(u'*Mu*u - 1) >= err
            c = u/norm(u);
            u = quadprog(2*G'*M*G, [], [], [], transpose(Mu*c), 1);
        end
    else
        G = grad(V,F);
        %M = speye(2*k);
        cellArea = repelem(power(doublearea(V, F)/2, 0.5), 2);
        M = spdiags(cellArea, 0, k*2, k*2);
        Mu = speye(n,n);
        %generate c here
        nullspace = null(transpose(Mu*u_prev));
        c = nullspace(:,1);
        c = c/norm(c);
        %generate linear constraint here
        A = zeros(max_iterations - 1, n);
        for i = 1:max_iterations
            if i == max_iterations
                A(i,1:n) = transpose(Mu*c);
            else
                A(i,1:n) = transpose(Mu*u_prev(:,i));
            end
        end
        b = [zeros(max_iterations - 1,1);1];
        u = quadprog(2*G'*M*G, [], [],[], A, b);
        while abs(u'*Mu*u - 1) >= err
            c = u/norm(u);
            for i = 1:max_iterations
                if i == max_iterations
                    A(i,1:n) = transpose(Mu*c);
                else
                    A(i,1:n) = transpose(Mu*u_prev(:,i));
                end
            end
            u = quadprog(2*G'*M*G, [],[], [], A, b);
        end
    end
end