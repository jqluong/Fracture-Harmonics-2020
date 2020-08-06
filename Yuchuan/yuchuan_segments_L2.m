%% Compute argmin_f ||f-g||^2 + ||f'||^2 using a line segment method

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
mu = 2.5;  %weight for the ||f-g||^2 term
lambda = 1;   %weight for the ||f'||^2 term
L = [eye(n-1) zeros(n-1,2) eye(n-1)];
M = mu*transpose(L)*L;
P = lambda*(1/h)*[-eye(n-1) zeros(n-1,2) eye(n-1)];
Q = transpose(P)*P;
D = [eye(n) -eye(n)];


%% Quadprog optimization
options = optimoptions(@quadprog, 'Algorithm', 'interior-point-convex', 'OptimalityTolerance', 1e-12);

f = quadprog(2*(M+Q), -2*transpose(g)*M, [], [], D, zeros(n,1),[],[],[], options); %with continuity constraint
%f = quadprog(2*(M+Q), -2*transpose(g)*M, [], [], [], [],[],[],[], options); %without continuity constraint


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





