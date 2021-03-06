
% for test: generate circle
s = 100;
theta = [0:2*pi()/s:2*pi()]';
x = cos(theta);
y = sin(theta);

[V, F] = triangle([x,y], 'Quality', 30, 'MaxArea', 0.001);

% for test: create square mesh
%x = [0:.1:1 1.5];
%y = x';
%[x,y] = meshgrid(x,y);
%x = x(:);
%y = y(:);

%F = delaunay(x,y);
%V = [x y];

% for test: create annulus
%[V,F,bo_in,bo_out] = annulus(100,2);


% set values for mu (see paper) and num (number of eigenmodes)
mu = 0;
num = 25;
U = zeros(length(V),num);      % matrix of eigenmodes, each column is an eigenmode.

% call function to find eigenmodes
U = eigenmodes_iterations(V,F,U, mu,num);

% plot eigenmode functions
filename = '/Users/lmattos/Documents/GitHub/Fracture-Harmonics-2020/Leticia/2D/vertex_func_eigenmodes.gif'; % Specify the output file name
plot_eigenmodes(V,F,U,filename);

function  R = eigenmodes_iterations(V,F,U, mu,num)

    CONVERG = 0.001;        % convegence criteria used to stop iterations

    L = -cotmatrix(V,F);    % laplacian matrix size #V by #V
    M = massmatrix(V,F);    % mass matrix, size #V by #V

    % find total of num eigenmodes loop
    for i = 1:num
        perform_iter(i);
    end  
    
    R = U; % set resulting matrix to U
    
    
    function perform_iter(i)
        c = rand(length(U),1);      % initialize the first vector to random unit vector
        c = c/(sqrt(c'*M*c));      % normalize the initial vector 
    
        % create inequality constraint
        Eq_0 = transpose(U(:,1:(i-1)))*M;     % matrix used for the orthogonal condition
        beq = [ zeros((i-1),1); 1];           % inequality constraint vector
    
        % find i-th eigenmode loop
        while sqrt(transpose(c-U(:,i))*M*(c-U(:,i)))  > CONVERG
            Eq_1 = transpose(c)*M;             % matrix used for the unit norm condition
            Aeq = [ Eq_0; Eq_1 ];   % append matrices for use in quadprog
            u = quadprog(L, mu*ones(1,length(V))*M, [],[],Aeq,beq);     % find solution
            U(:,i) = c;     % store solution into matrix of eigenmodes
            c = u/(sqrt(u'*M*u));  % normalize solution
        end % end of find i-th eigenmode loop
    
    end % end of perform_iter function
    
end % end of eigenmodes_iterations function

function plot_eigenmodes(V, F, U, filename)
    
    h = figure; 
    [~,nImages] = size(U);
    axis tight manual;
    delay = 0.25;
    
    for idx = 1:nImages
        % capture plot sequence of images
        tsurf(F,[V U(:,idx)],fphong, falpha(1,0));
        view(0,90);
        axis off;
        axis equal;
        title(['Eigenmodes of L using quadprog,  n = ' num2str( idx) ])
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