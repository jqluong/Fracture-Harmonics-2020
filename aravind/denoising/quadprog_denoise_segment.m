%% problem
% argmin_f |f-g|_2^2 + alpha*|f'|_2 + beta*|Df|_1
alpha = 1e6;
beta = 10;

% construct noisy function
n = 1000; % number of segments, (divisible by 2 for easy function construct.)
nf = 100; % frequency of noise

x = 1/n*[0:n-1, 1:n]';
gi = kron([1;0;1;0], ones(n/2,1)); % no noise
g = gi;
for i=2:n
    g(i) = gi(i) + 0.1*sin(nf*i/n);
    g(i+n-1) = g(i);
end

f = quadprog_denoise(g, alpha, beta);

segment_metric(x,f,gi)

% plot results
hold on
segment_plot(x,gi,'c')
segment_plot(x,g,'b')
segment_plot(x,f, 'r')
segment_plot_legend({'initial signal','noisy signal','output signal'}, ['c','b','r']);
title(['Denoising with quadprog; \alpha = ', num2str(alpha), ', \beta = ', num2str(beta)]);
hold off

%% quadprog function 
function f = quadprog_denoise(g, alpha, beta)
    n = length(g)/2;
    h = 1/n;
    
    % base matrices
    G = [segment_derivative(n), zeros(n)];
    D = [beta*segment_discontinuity(n); zeros(1, 2*n)]; 
    %M = [speye(2*n), zeros(2*n, n); zeros(n, 3*n)];
    L = segment_midpoint(n);
    L = 2*(L'*L);
    M = [L, zeros(2*n, n); zeros(n, 3*n)];

    % quadratic
    H = alpha*(G'*G) + M;

    % linear
    %v = [-2*g', ones(1,n)];
    v = [-2*g'*L, ones(1,n)];

    % constraint
    C = [-D, -speye(n); D, -speye(n)];
    z = zeros(size(C,1),1);

    % Minimization
    opts = optimoptions(@quadprog, 'Algorithm', 'interior-point-convex', 'OptimalityTolerance', 1e-6, 'Display', 'off');
    ft = quadprog(2*H,v,C,z,[],[],[],[],[],opts);
    f = ft(1:2*n);
end