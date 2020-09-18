%% Compute eigenmodes of Laplacian on mesh for vertex functions

%% Create mesh
%f = 20;
%x = 0:1/f:1;
%y = 0:1/f:2;

%[X,Y] = meshgrid(x,y);
%points = [X(:), Y(:)];

%V = points;
%F = delaunay(points);

k = 512; % boundary spacing: 2pi/n
th = [0:2*pi()/k:2*pi()]';
x = cos(th);
y = sin(th);

[V, F] = triangle([x,y], 'Quality', 30, 'MaxArea', 0.001);

%% Initialization 
n = 25; % number of eigenmodes to find
U = zeros(length(V), n);

%% Find eigenmodes
U = laplacian_eigenmodes_iterative(V, F, U, n);

%% Plot eigenmodes
eigenmodes_gif(V, F, U, 'C:\Users\aravi\Desktop\laplacianeigenmodes.gif');

%% Iterative method function
function Z = laplacian_eigenmodes_iterative(V, F, U, n)
    % V, F are mesh
    % U is the matrix to store eigenmodes
    % n is the number of iterations, matches #columns of U
    
    TOL = 0.0001;
    opts = optimoptions(@quadprog, 'Algorithm', 'interior-point-convex', 'OptimalityTolerance', 1e-6, 'Display', 'off');

    G = grad(V, F);
    A = sparse(1:size(G,1), 1:size(G,1), 1/6*kron([1;1],doublearea(V,F)));
    MG = A.^0.5 * G;
    GMG = MG'*MG; % quadratic optim term
    M = massmatrix(V,F);
    
    for i = 1:n % find n eigenmodes
        disp(i);
        iterative_computation(i);
    end
    
    Z = U;
    
    function iterative_computation(i)
        c = rand(length(U), 1);
        c = c / norm(c);
        Ai = (U(:,1:(i-1))')*M; % linear equality constraint - orthgonal to previous solutions
        beq = [zeros((i-1), 1); 1]; % linear equality constraint

        %while abs(1 - c'*M*U(:,i)) > TOL % this does not work :\
        while norm(c - U(:,i)) > TOL
            Aeq = [Ai; c'*M]; % lienar equality constraint - append "unit norm"
            u = quadprog(GMG, [], [], [], Aeq, beq, [], [], [], opts);
            U(:,i) = c;
            c = u/norm(u);
        end
    end
end

function eigenmodes_gif(V, F, U, fn)
    m = max(max(abs(U)));
    n = size(U, 2);
    h = figure;  
    set(gcf,'Color','w');
    axis tight manual;
    delay = 0.5; % gif delay time
    
    for i = 1:n
        % draw plot
        tsurf(F,[V U(:,i)],fphong, falpha(1,0));
        set_display(i);
        drawnow;
        % capture plot as an image 
        frame = getframe(h);
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256);
        % write to gif file
        if i == 1 
            imwrite(imind,cm,fn,'gif', 'Loopcount',inf, 'DelayTime', delay); 
        else 
            imwrite(imind,cm,fn,'gif','WriteMode','append', 'DelayTime',delay); 
        end 
    end
    
    function set_display(i)
        zlim([-m m]);
        caxis([-m m]);
        colormap(cbrewer('RdYlBu',40));
        colorbar();
        title(strcat('Laplacian Eigenmode', {' '}, num2str(i)));        
        grid off
        axis off
        axis equal
        set(gcf,'Color','w');
    end
end