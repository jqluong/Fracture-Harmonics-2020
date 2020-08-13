%% Compute argmin_f mu*||f-g||^2 + lambda*||f'||^2 + eta*||Df||^2 using a line segment method

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
mu = 1;  %weight for the ||f-g||^2 term
lambda = 1;   %weight for the ||f'||^2 term
eta = 1000; %weight for the continuity term
L = sparse([speye(n-1) zeros(n-1,2) speye(n-1)]);
M = mu*transpose(L)*L;
P = sparse((1/h)*[-speye(n-1) zeros(n-1,2) speye(n-1)]);
Q = lambda*transpose(P)*P;
D = sparse([speye(n) -speye(n)]);
DTD = eta*transpose(D)*D;


%% Quadprog optimization
options = optimoptions(@quadprog, 'Algorithm', 'interior-point-convex', 'OptimalityTolerance', 1e-8);

%f = quadprog(2*(M+Q), -2*transpose(g)*M, [], [], D, zeros(n,1),[],[],[], options); %with continuity constraint
%f = quadprog(2*(M+Q), -2*transpose(g)*M, [], [], [], [],[],[],[], options); %without continuity constraint
f = quadprog(2*(M+Q+DTD), -2*transpose(g)*M, [], [], [], [],[],[],[], options); %with continuity term


%% Plotting results

figure(1) %plot line segments

hold on

for i = 1:n-1
    plot([t(i) t(i+1)], [f(i) f(i+1+n)],'LineWidth',1,'color','blue');
end

for i = 1:n-1
     plot([t(i) t(i+1)], [g(i) g(i+1+n)],'LineWidth',1,'color',[0.8500, 0.3250, 0.0980]);
end


hold off





