function  eigenmodes_quadprog(V,F,mu,num)

CONVERG = 0.001;        % convegence criteria used to stop iterations

[n,~] = size(V);

[m,~] = size(F);

L = -cotmatrix(V,F);    % laplacian matrix size #V by #V
M = massmatrix(V,F);    % mass matrix, size #V by #V
U = zeros(n,num);      % matrix of eigenmodes, each column is an eigenmode.
[a,~] = size(U);

% find total of num eigenmodes loop
for i = 1:num
    j = i-1;
    
    c = rand(a,1);      % initialize the first vector to random unit vector
    c = c/norm(c);      % normalize the initial vector
    
    % create inequality constraint
    Eq_0 = transpose(U(:,1:j))*M;     % matrix used for the orthogonal condition
    Eq_1 = transpose(c)*M;             % matrix used for the unit norm condition
    Aeq = [ Eq_0; Eq_1 ];   % append matrices for use in quadprog
    beq = [ zeros(j,1); 1]; % inequality constraint vector
    
    diff = norm(c-U(:,i));
    
    % find i-th eigenmode loop
    while diff > CONVERG
        u = quadprog(L, mu*ones(1,n)*M, [],[],Aeq,beq);     % find solution
        c = u/norm(u);  % normalize solution
        U(:,i) = c;     % store solution into matrix of eigenmodes
    end % end of find i-th eigenmode loop
    
end % end of find eigenmodes loop

% plot figures

nImages = length(U);

fig = figure;
for idx = 1:nImages
    tsurf(F,[V U(:,idx)],fphong, falpha(1,0));
    axis equal;
    set(gcf,'Color','w'); 
    delay = 5;
    title(['eigenmodes using quadprog,  n = ' num2str( n(idx)) ])
    drawnow
    frame = getframe(fig);
    im{idx} = frame2im(frame);
end
close;

end