%Finding eigenvalues for the Laplacian by iteratively minimizng Dirichlet
%Energy.

%Borrowing functions from this folder
mkdir('vertex_functions')
addpath('vertex_functions')

%Square Mesh
x = [0:.1:1 1.5];
y = x';
[x,y] = meshgrid(x,y);
x = x(:);
y = y(:);

%Circular Mesh
%s = 100;
%theta = [0:2*pi()/s:2*pi()]';
%x = cos(theta);
%y = sin(theta);

%Generate Mesh
F = delaunay(x,y);
V = [x y];

%Set number of iterations
max_iterations = 3;
[m,~] = size(V);
%Store u values into this matrix where u(:,i) is the ith eigenfunction
u_prev = zeros(m, max_iterations);
for i = 1:max_iterations
    u = laplace_eq_2D_quadprog_iterations(V, F, u_prev, i, max_iterations);
    u_prev(:,i) = u;
end
%how to plot triangle function (reminder, for me):
%plot_surf(V,F,u)


function u = laplace_eq_2D_quadprog_iterations(V, F, u_prev, iterations, max_iterations)
    % solve laplace's equation via quadprog.  Obtain eigenfunctions by
    % adding the constraints that u^T*M*u = 1 and u_i^T*M*u_j = 0 for all j
    % less than i.
    
    [n,~] = size(V);
    [k, ~] = size(F);
    err = 10^-5;
    if iterations == 1
        %Setting up matrices
        %Define gradient matrix
        G = grad(V,F);
        %Define M matrix for G^TMG
        %M = speye(2*k);
        cellArea = repelem(power(doublearea(V, F)/2, 0.5), 2);
        M = spdiags(cellArea, 0, k*2, k*2);
        %dblA = doublearea(V,F);
        %M = repdiag(diag(sparse(dblA)/2),size(V,2));
        %Define a weighting matrix for u^T*M*u (weighted dot product)
        Mu = speye(n,n);
        
        %Create a random unit vector to simulate u^T*M*u = 1
        c = rand(n,1);
        c = c/norm(c);
        
        %Keep solving until u^T*M*u is close enough to 1
        u = quadprog(2*G'*M*G, [], [], [], transpose(Mu*c), 1);
        while abs((u'*Mu*u) - 1) >= err
            %If the "l2 norm" of u is not close enough to 1, set c to
            %previous u and try again
            c = u/norm(u);
            u = quadprog(2*G'*M*G, [], [], [], transpose(Mu*c), 1);
        end
    else
        %This is the case where you need to enforce orthogonality
        %constraints on u
        
        G = grad(V,F);
        %M = speye(2*k);
        cellArea = repelem(power(doublearea(V, F)/2, 0.5), 2);
        M = spdiags(cellArea, 0, k*2, k*2);
        Mu = speye(n,n);
        %dblA = doublearea(V,F);
        %M = repdiag(diag(sparse(dblA)/2),size(V,2));
        
        %Generate a random unit vector c that is orthongal to the previous
        %u's
        nullspace = null(transpose(Mu*u_prev));
        c = nullspace(:,1);
        c = c/norm(c);
        %Generate linear constraint here
        A = zeros(max_iterations, n);
        for i = 1:max_iterations
            if i == max_iterations
                %Enforce l2 constraint
                A(i,:) = transpose(Mu*c);
            else
                %Enforce orthogonality constraint
                A(i,:) = transpose(Mu*u_prev(:,i));
            end
        end
        b = [zeros(max_iterations - 1,1);1];
        %Keep solving until u^TMu is close enough to 1
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