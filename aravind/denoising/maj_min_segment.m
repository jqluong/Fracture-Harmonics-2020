%% problem
n=100; % number of segments
t=0:1/n:1;

x = repelem(t,[1, 2*ones(1,n-1), 1]); % segment positions

% initial signal for reference
yi1 = ones(1,n/2+1);
yi2 = zeros(1, n/2+1);
yi = [repelem(yi1, [1, 2*ones(1, length(yi1)-2), 1]), repelem(yi2, [1, 2*ones(1, length(yi2)-2), 1])]';

% construct randomized signal
y1 = ones(1,n/2+1) + 0.2*(rand(1,n/2+1) - 0.5);
y2 = zeros(1, n/2+1) + 0.2*(rand(1,n/2+1) - 0.5);
y = [repelem(y1, [1, 2*ones(1, length(y1)-2), 1]), repelem(y2, [1, 2*ones(1, length(y2)-2), 1])]';

f = mm_denoise(y);

NRMSE = sqrt(mean((yi - f).^2))/range(yi);
fprintf('normalized root mean square error for l1 signal recovery: %f \n', NRMSE);

plot_signals(x, y, f);

%% Majorization-Minimization Algorithm
function x = mm_denoise(y)
    % implement majorization-minimization algorithm to find x that minimizes
    % 0.5 ||x - y||^2 + \lambda |x'| 
    %
    % y is input noisy signal
    %
    % formulation of the problem is described in: 
    % http://eeweb.poly.edu/iselesni/lecture_notes/TVDmm/TVDmm.pdf
    
    % behavioral variables
    MAX_ITER = 1000;
    TOL = 0.001;
    ERR = 0.001;
    
    % algorithm variables
    n = length(y)/2;
    D = discrete_derivative_segment(n,1);
    lambda = 0.5;
    DDT = D*D';
    x = 0.5*ones(length(y),1);
    z = x;
    Dx = D*x + ERR;
    Dy = D*y;
    
    % majorization-minimization
    for k = 1:MAX_ITER
        z = x;
        A = sparse(1/lambda*(diag(abs(Dx)) + ERR)) + DDT;
        x = y - D'*(A\Dy);
        
        if sqrt(mean((x - z).^2))/range(x) < TOL
            break;
        end
        
        Dx = D*x;
    end
end

function plot_signals(x, s1, s2)
    x1 = x(1:2:end); % line segment position 1
    x2 = x(2:2:end); % line segment position 2
    
    s1a = s1(1:2:end)';
    s1b = s1(2:2:end)';
    
    s2a = s2(1:2:end)';
    s2b = s2(2:2:end)';
    
    hold on
    plot([x1;x2],[s1a;s1b], 'Color', 'b');
    plot([x1;x2],[s2a;s2b], 'Color', 'r');
    %legend('Noisy Signal', 'Recovered Signal')
    hold off
end