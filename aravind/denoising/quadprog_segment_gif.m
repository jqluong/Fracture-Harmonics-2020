%% problem
% argmin_f |f-g|_2^2 + alpha*|f'|_2 + beta*|Df|_1
alpha = 0;
beta = 0;

fig = figure;
%axis tight manual
filename = 'QuadProgAnimated.gif';

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
m = segment_metric(x,f,gi);

        hold on
        segment_plot(x,gi,'c')
        segment_plot(x,g,'b')
        segment_plot(x,f, 'r')
        segment_plot_legend({'initial signal','noisy signal','output signal'}, ['c','b','r']);
        title(['\alpha = ', num2str(alpha), ', \beta = ', num2str(beta), ', metric = ', num2str(m)]);
        hold off
        drawnow;
        frame = getframe(fig); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);

for alpha = [10, 500, 1000, 5000, 10000, 100000, 1000000]
    for beta = [1, 5, 10, 20, 100]        
        f = quadprog_denoise(g, alpha, beta);
        m = segment_metric(x,f,gi);
        
        % plot results
        clf(fig,'reset')
        hold on
        segment_plot(x,gi,'c')
        segment_plot(x,g,'b')
        segment_plot(x,f, 'r')
        segment_plot_legend({'initial signal','noisy signal','output signal'}, ['c','b','r']);
        title(['\alpha = ', num2str(alpha), ', \beta = ', num2str(beta), ', metric = ', num2str(m)]);
        hold off
        drawnow;
        
        frame = getframe(fig); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    end
end

%% quadprog function 
function f = quadprog_denoise(g, alpha, beta)
    n = length(g)/2;
    h = 1/n;
    
    % base matrices
    G = [segment_derivative(n), zeros(n)];
    D = [beta*segment_discontinuity(n); zeros(1, 2*n)];
    I = [speye(2*n), zeros(2*n, n); zeros(n, 3*n)];

    % quadratic
    H = alpha*(G'*G) + I;

    % linear
    v = [-2*g', ones(1,n)];

    % constraint
    C = [-D, -speye(n); D, -speye(n)];
    z = zeros(size(C,1),1);

    % Minimization
    opts = optimoptions(@quadprog, 'Algorithm', 'interior-point-convex', 'OptimalityTolerance', 1e-6, 'Display', 'off');
    ft = quadprog(2*H,v,C,z,[],[],[],[],[],opts);
    f = ft(1:2*n);
end
