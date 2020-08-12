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

f = mm_denoise(s);

NRMSE = sqrt(mean((si - f).^2))/range(si);
fprintf('normalized root mean square error for l1 signal recovery: %f \n', NRMSE);

hold on
plot(0:1/(n-1):1, s);
plot(0:1/(n-1):1, f);
legend('Noisy Signal', 'Recovered Signal')
hold off

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
    
    % algorithm variables
    n = length(y);
    D = sparse(discrete_derivative(n,1));
    lambda = 0.5;
    DDT = D*D';
    x = 0.5*ones(n,1);
    z = 0.5*ones(n,1);
    Dx = D*x;
    Dy = D*y;
    
    % majorization-minimization
    for k = 1:MAX_ITER
        z = x;
        A = sparse(1/lambda*diag(abs(Dx))) + DDT;
        x = y - D'*(A\Dy);
        
        if sqrt(mean((x - z).^2))/range(x) < TOL
            break;
        end
        
        Dx = D*x;
    end
end