%% Compute argmin_f mu*||f-g||^2 + lambda*||f'|| + eta*||Df|| using a line segment method
clear

%% Input data
a = 0; %left endpoint
b = 10; %right endpoint
t = (a:0.05:b)';   %partition of [a,b]
s = 10*square((pi/5)*t);         %square wave 
%s = 10*sawtooth(t);             %sawtooth wave
noisy = awgn(s,10,'measured');
g = [noisy;noisy];                %concatenated 2 copies of the same noisy signal
n = length(t);                    %number of line segments
h = (b-a)/n;


%% Build matrices

%weights
mu = 1;  %weight for the ||f-g||^2 term
lambda = 1;   %weight for the ||f'||^2 term
eta = 100; %weight for the continuity term

%matrices for mu*||f-g||^2 term
L = sparse([speye(n-1) zeros(n-1,2) speye(n-1)]);
M = sparse(mu*transpose(L)*L);

%matrices for lambda*||f'|| term
P = sparse(lambda*(1/h)*[-speye(n-1) zeros(n-1,2) speye(n-1)]);
Q = lambda*transpose(P)*P;

%matrices for eta*||Df|| term
D = sparse(eta*[speye(n) -speye(n)]);

%matrix for quadratic term
H = sparse(3*n,3*n);
H(1:2*n,1:2*n) = M+Q;
v = ones(3*n,1);
v(1:2*n) = -2*M*g;

%inequality constraint matrix
ineq = sparse([D -1*speye(n) ; -D -1*speye(n)]);


%% Minimization
options = optimoptions(@quadprog, 'Algorithm', 'interior-point-convex', 'OptimalityTolerance', 1e-8);
ft = quadprog(2*H,v,ineq,zeros(2*n,1));
f = ft(1:2*n);

%% Plotting results

figure(1) %plot line segments
hold on

for i = 1:n-1
    plot([t(i) t(i+1)], [f(i) f(i+1+n)],'color','blue');
end

for i = 1:n-1
     plot([t(i) t(i+1)], [g(i) g(i+1+n)],'color',[0.8500, 0.3250, 0.0980]);
end

hold off

