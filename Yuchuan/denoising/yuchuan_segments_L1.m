%% Compute argmin_f ||f-g||^2 + ||f'|| using a line segment method

clear


%% Input data
a = 0; %left endpoint
b = 10; %right endpoint
t = (a:0.05:b)';
sq = 10*square((pi/5)*t);         %square wave of period = 10            
noisy = awgn(sq,10,'measured');
g = [noisy;noisy];                %concatenated 2 copies of the same noisy signal
n = length(t);                    %number of line segments
h = (b-a)/n;

%% Build matrices

L = [eye(n-1) zeros(n-1,2) eye(n-1)];
M = transpose(L)*L;
P = (1/h)*[-eye(n-1) zeros(n-1,2) eye(n-1)];
D = [eye(n) -eye(n)];






%% TV Denoising
% Run TV denoising algorithm (MM algorithm)

lam = 2;                         % lam: regularization parameter
Nit = 500;                          % Nit: number of iterations
f = tvdmm_segments(g, lam, Nit, h);   % Run MM TV denoising algorithm






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


