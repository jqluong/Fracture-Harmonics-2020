function face_min_dirichlet_L1(V,F, num, filename)

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

% plot eigenmode functions
animated_eigenmodes(Vd,Fd,k, Y, filename)

function  R = eigenmodes_iterations(Vd,Fd,D,m,k,Y,num)

    CONVERG = 0.001;          % convegence criteria used to stop iterations
    beta = 1e-7;

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


function animated_eigenmodes(V,F,k, Y, filename)
    
    h = figure; 
    [~,nImages] = size(Y);
    axis tight manual;
    delay = 1;
    
    for idx = 1:nImages
        % capture plot sequence of images
        face_plotting(V,F,Y(1:3*k,idx));
        zlim([-0.5 0.5]);
        caxis([-0.5 0.5]);
        axis off;
        colormap(cbrewer('RdYlBu',40));
        colorbar();
        view(345,45);
        title(['Iterative min_{u} 1/2*u^T*L*u + 1^T*|Du|,  n = ' num2str( idx) ])
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


end