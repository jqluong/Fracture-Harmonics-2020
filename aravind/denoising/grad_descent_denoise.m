%% problem
freq=1000;
t=0:1/freq:1;

si=zeros(length(t)-1,1);
for i = 1:(length(si)/2)
    si(i) = 1;
end
sr = 0.2*(rand([length(si),1])-0.5);
s=si+sr;
n = length(s);

f = gradient_descent_denoise(s);

NRMSE = sqrt(mean((si - f).^2))/range(si);
fprintf('normalized root mean square error for l1 signal recovery: %f \n', NRMSE);

hold on
plot(0:1/(n-1):1, s);
plot(0:1/(n-1):1, f);
legend('Noisy Signal', 'Recovered Signal')
hold off

%% gradient descent
function x = gradient_descent_denoise(y)
    % implement gradient descent algorithm to find x that minimizes
    % 0.5 ||x - y||^2 + \lambda |x'| 
    % y is input noisy signal
    
    % behavioral variables
    %tstart=tic;
    MAX_ITER = 100;
    TOL = 0.01;
    
    % algorithm variables
    n = length(y);
    D = discrete_derivative(n,1); % fails dramatically when h is correct?
    dx = zeros(n,1); % vector for discrete gradient
    ds = zeros(n,1); % initialize
    du = zeros(n,1); % initialize
    t = 1;
    lambda = 0.5;
    x = 0.5*ones(n,1); 
    a = 0.2; % backtracking line search parameter a in (0,0.5)
    b = 0.3; % backtracking line search parameter b in (0,1)
    
    % gradient descent
    for k = 1:MAX_ITER
        
        % first, find dx = - del_x (f)
        dx = -del_f(x);
        ndx = norm(dx);
        
        % check if finished
        if ndx < TOL
            break; 
        end
        
        % now, use backtracking line search to find step size t
        t = 1;
        while objective(x + t*dx) > objective(x) - a*t*norm(dx)^2
            t = b*t;
        end
        
        %update
        x = x + t*dx;
    end
    %toc(tstart);
    
    function m = objective(s)
        m = 0.5*norm(s-y)^2 + lambda*norm(D*s,1);
    end
    
    function dx = del_f(s)
        % del_x f = del_x (0.5 || x - y ||_2^2 + \lambda || x' ||_1)
        
        ds = D*s; % x'
        % when normalizing, make sure no divide by 0
        du = D*(ds ./ sqrt((0.001^2 + ds.^2)));
        dx = (s - y) - lambda*du;     
    end
end