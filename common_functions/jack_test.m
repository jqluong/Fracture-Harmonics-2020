%jack's testing script, it might get moved to my folder later
%works with normal old face functions
mkdir('vertex_functions')
addpath('vertex_functions')
mkdir('face_functions')
addpath('face_functions')
mkdir('mesh_functions')
addpath('mesh_functions')
%Testing circle meshes
%[V,F] = mesh_circle(10);
%M = face_build_discontinuity_mass(V,F);
%L = face_build_discontinuity_laplacian(V,F);
%[m,n] = size(F);
%f = ones(3 * m, 1);
%M * f
%L * f


x = 0:.1:1;
y = x';
[x,y] = meshgrid(x,y);
x = x(:);
y = y(:);
F = delaunay(x,y);
V = [x y];
max_iterations = 2;
[m,~] = size(V);
u = laplace_disc(V,F)

%get eigenfunctions
%u_prev = zeros(m*(max_iterations - 1), 1);
%for i = 1:max_iterations
%    u = laplace_eq_2D_quadprog_iterations(V, F, u_prev, i);
%    u_prev(i*(m-1) + 1: m*i) = u;
%end
%u1 = laplace_eq_2D_seg_quadprog(V, F);
%face_plotting(V,F,u1)
%U = iterative(V,F,max_iterations);

function Y = iterative(V,F, num)

% solves min_{u} 1/2*u^T*L*u + 1^T*|Du|, where u is a face function, using
% iterative method and quadprog
%
% now it outputs results into GIF file
%
% V: continuous mesh vertex matrix, size #|V| by 2
% F: continuous mesh face matrix, size #|F| by 3
% num: number of eigenmodes desired
% filename: string, specify file directory to write GIF file
%
% u: discontinuous mesh solution vector, size #3|F| by 1


% transform mesh matrices to discontinuous mesh space
% size |Fd|=|Fc|
% size |Vd| is variable
[Vd,Fd] = discontinuous_reshape(V,F);

% declare some dimensions used
E = edges(F);        % edges matrix
[m,~] = size(E);     % needed for dimension of Du
[k,~] = size(F);     % needed for length u in R^3*k

Y = zeros(3*k,num); % matrix of eigenmodes 3*|F| by num, each column is an eigenmode.

% generate discontinuity matrix, size 2*|E| by 3*|F|
D = discontinuity([V zeros(length(V),1)], F);

% call function to find eigenmodes
Y = eigenmodes_iterations(Vd,Fd,D,m,k,Y,num);


function  R = eigenmodes_iterations(Vd,Fd,D,m,k,Y,num)

    CONVERG = 0.001;          % convegence criteria used to stop iterations
    beta = 1;

    L = -cotmatrix(Vd,Fd);    % disc. laplacian matrix size #3|F| by #3|F|
    M = massmatrix(Vd,Fd);    % mass matrix, size #3|F| by #3|F|
    
    % construct block matrix using Laplacian to act on y=[u;t]
    % recall first min term is 1/2 transpose(y)*L_tilde*y
    O_r = sparse(3*k,2*m);
    O_b = sparse(2*m,3*k+2*m);
    L_tilde = [L O_r; O_b];
    
    % construct block matrix using discontinuity matrix to act on y=[u;t]
    I = speye(2*m);
    D_tilde = beta * [-D -I; D -I];

    % construct vector f used for L1 norm
    % recall second min term is transpose(f)*y
    u_placeholder = zeros(3*k,1);
    t_placeholder = ones(2*m,1);
    f = [ u_placeholder; t_placeholder ];

    % find total of num eigenmodes loop
    for i = 1:num
        perform_iter(i);
        % display progress
        fprintf(['at ' num2str(i) '-th iteration'])
    end  
    
    R = Y; % set resulting matrix to U

    
    function perform_iter(i)
        c = rand(3*k,1);        % initialize the first vector to random unit vector
        c = c/(sqrt(c'*M*c));   % normalize the initial vector 
    
        % create equality constraint
        Eq_0 = transpose(Y(:,1:(i-1)))*M;     % matrix used for the orthogonal condition
        beq = [ zeros((i-1),1); 1];           % inequality constraint vector
    
        % find i-th eigenmode loop
        while sqrt(transpose(c-Y(:,i))*M*(c-Y(:,i)))  > CONVERG
            % matrix used for the unit norm condition
            Eq_1 = transpose(c)*M;  
            % append matrices for equality constraint
            Aeq = [ Eq_0 sparse(length(Eq_0(:,1)),2*m); Eq_1 sparse(length(Eq_1(:,1)),2*m)];
            % find solution using quadprog
            u = quadprog(L_tilde, f,D_tilde ,zeros(2*2*m,1),Aeq,beq);     
            Y(:,i) = c;     % store solution into matrix of eigenmodes
            c = u(1:3*k)/(sqrt(u(1:3*k)'*M*u(1:3*k)));  % normalize solution
            
        end % end of find i-th eigenmode loop
    
    end % end of perform_iter function
    
end % end of eigenmodes_iterations function
end

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
    O_right = zeros(3*k,m); 
    O_bottom = zeros(m,3*k+m);
    G_tilde = [L O_right; O_bottom ];
    
    % generate discontinuity matrix
    V = [V zeros(length(V))];
    D = discontinuity(V,F); %D = 2|E|-by-3|F| sparse matrix , E = #edges
    D_sum = face_edge_sum(F); %D_sum = 2|E| by |E| matrix
    V = V(:,1:2);

    
    % generate vector f used in the constraint
    u_placeholder = zeros(3*k,1);
    t_placeholder = ones(m,1);
    f = [ u_placeholder; t_placeholder ];
    
    % generate |E|x|E| identity matrix used in the constraint block matrix
    I = speye(m);
    
    % generate inequality constraint matrix and vector
    A = [D_sum*D -I; -D_sum*D -I];
    b = zeros(2*m,1);
    
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
    
    equalityConstraint = sparse(BP2new,ones(length(BP2new),1),ones(length(BP2new),1), 3*k + m, 1);
    equalityMatrix = sparse([BP1new; BP2new], [BP1new; BP2new], ones(length(BP1new) + length(BP2new), 1), 3*k + m, 3*k + m);
    
    y = quadprog(2*G_tilde,f,A,b,equalityMatrix,equalityConstraint);
    
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
%Comparing energy of quadprog solution v.s. "step function"
norm(L*v,2)^2 + norm(D_sum*D*v,1)
norm(L*u,2)^2 + norm(D_sum*D*u,1)
end