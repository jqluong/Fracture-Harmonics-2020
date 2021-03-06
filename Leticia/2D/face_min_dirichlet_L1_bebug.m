
function Y = face_min_dirichlet_L1_bebug(V,F, num)

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


% generate discontinuity matrix, size 2*|E| by 3*|F|
D = discontinuity([V sparse(length(V),1)], F);
%D_sum = face_edge_sum(F);
%D_new = D_sum * D;
V = V(:,1:2);

% transform mesh matrices to discontinuous mesh space
% size |Fd|=|Fc|
% size |Vd| is variable
[Vd,Fd] = discontinuous_reshape(V,F);

% declare some dimensions used
E = edges(F);        % edges matrix
[m,~] = size(E);     % needed for dimension of Du
[k,~] = size(F);     % needed for length u in R^3*k

Y = sparse(3*k,num); % matrix of eigenmodes 3*|F| by num, each column is an eigenmode.z



% call function to find eigenmodes
Y = eigenmodes_iterations(Vd,Fd,D,m,k,Y,num);


filename = '/Users/lmattos/Documents/GitHub/Fracture-Harmonics-2020/Leticia/2D/2D_pictures/animated_eigen_L1_beta_1e2_disregard.gif';

% plot eigenmode functions
animated_eigenmodes(V,F, full(Y), filename);



function  R = eigenmodes_iterations(Vd,Fd,D,m,k,Y,num)

    CONVERG = 0.001;          % convegence criteria used to stop iterations
    beta = 1e-3;

    L = -cotmatrix(Vd,Fd);    % disc. laplacian matrix size #3|F| by #3|F|
    M = massmatrix(Vd,Fd);    % mass matrix, size #3|F| by #3|F|
    
    % construct block matrix using Laplacian to act on y=[u;t]
    % recall first min term is 1/2 transpose(y)*L_tilde*y
    O_r = sparse(3*k,2*m);
    O_b = sparse(2*m,3*k+2*m);
    L_tilde = [L O_r; O_b];
    
    % construct block matrix using discontinuity matrix to act on y=[u;t]
    I = speye(2*m);
    D_tilde = [-D -I; D -I];

    % construct vector f used for L1 norm
    % recall second min term is transpose(f)*y
    u_placeholder = sparse(3*k,1);
    t_placeholder = ones(2*m,1);
    f = [ u_placeholder; t_placeholder ];

    % find total of num eigenmodes loop
    for i = 1:num
        perform_iter(i);
        % display progress
        fprintf(['at ' num2str(i) '-th iteration \n'])
    end  
    
    R = Y; % set resulting matrix to U

    
    function perform_iter(i)
        options =  optimset('Display','off');
        c = rand(3*k,1);        % initialize the first vector to random unit vector
        c = c/(sqrt(c'*M*c));   % normalize the initial vector 
    
        % create equality constraint
        Eq_0 = transpose(Y(:,1:(i-1)))*M;     % matrix used for the orthogonal condition
        beq = [ sparse((i-1),1); 1];           % inequality constraint vector
    
        % find i-th eigenmode loop
        while sqrt(transpose(c-Y(:,i))*M*(c-Y(:,i)))  > CONVERG
            % matrix used for the unit norm condition
            Eq_1 = transpose(c)*M;  
            % append matrices for equality constraint
            Aeq = [ Eq_0 sparse(length(Eq_0(:,1)),2*m); Eq_1 sparse(length(Eq_1(:,1)),2*m)];
            % find solution using quadprog
            u = quadprog(2*L_tilde, beta*f,D_tilde ,sparse(4*m,1),Aeq,beq,[],[],[],options);     
            Y(:,i) = c;     % store solution into matrix of eigenmodes
            c = u(1:3*k)/(sqrt(u(1:3*k)'*M*u(1:3*k)));  % normalize solution
            
        end % end of find i-th eigenmode loop
    
    end % end of perform_iter function
    
end % end of eigenmodes_iterations function

end

function step = create_step(V,F)
    [m,~] = size(F);
    step = -1/10*ones(3 * m,1);
    for i = 1:m
        vertices = V(F(i,:));
        if vertices(1) >= 5 && vertices(2) >= 5 && vertices(3) >= 5
            step(3*(i-1) + 1) = 1/10;
            step(3*(i-1) + 2) = 1/10;
            step(3*(i-1) + 3) = 1/10;
        end
    end
end

function [dirichlet, disE] = compute_energy(V,F,u)
    L = -face_build_discontinuity_laplacian(V,F);
    [m,~] = size(V);
    V = [V zeros(m,1)];
    D = discontinuity(V,F);
    D_sum = face_edge_sum(F);
    D_new = D_sum * D;
    %V = V(:,1:2);
    dirichlet = (norm(L*u,2)^2);
    disE = norm(D_new*u,1);
end

function animated_eigenmodes(V,F, Y, filename)
    
    h = figure; 
    [~,nImages] = size(Y);
    axis tight manual;
    delay = 1;
    
    for idx = 1:nImages
        % capture plot sequence of images
        face_plotting(V,F,Y(:,idx));
        zlim([-0.5 0.5]);
        caxis([-0.5 0.5]);
        axis off;
        colormap(cbrewer('RdYlBu',40));
        colorbar();
        view(-115,45);
        title(['Iterative min_{u} alpha*1/2*u^T*L*u + beta*1^T*|Du|,  n = ' num2str( idx) ' alpha = 2, beta = 1e-3'])
        drawnow;
        frame = getframe(h);
        im = frame2im(frame); 
        [A,map] = rgb2ind(im,256);
        % save images into a GIF file
        if idx == 1 
            imwrite(A,map,filename,'gif', 'Loopcount',inf, 'DelayTime', delay); 
        else 
            imwrite(A,map,filename,'gif','WriteMode','append', 'DelayTime',delay); 
        end 
    end
    
    


end

