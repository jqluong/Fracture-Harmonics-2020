function min_quad_prog(numsteps)

% example:
%
% 0---x---x---x---x---1
%
% n = 5
%
% number of function values is 6

% number of divisions on interval [0,1]
n = numsteps;

% create noisy input 
% first use standard square with desired amplitude and period
x = 0:2*pi/n:2*pi;
y = 0.5+ 0.25*square(x);
% then add random noise
g0 = y + (0.1) * randn(1,n+1);

% example:
%
% g0 = 0.81524      0.78271      0.85826      0.35061      0.18491      0.77571
%
% g  = 0.81524      0.78271      0.85826      0.35061      0.18491      0.78271      0.85826      0.35061      0.18491      0.77571


% keep first n values
g = zeros(2*n,1);
for i = 1:n
    g(i) = g0(i);
end
% concatenate additional values
for i = n+1:2*n
    g(i) = g0(i-n+1);
end

% gradient matrix
% G(i,i) = 1
% G(i,i+n) = -1
%
% example:
%
% i = [ 1 2 3 ... 5  1  2  3 ...  5 ]
%
% j = [ 1 2 3 ... 5  6  7  8 ... 10 ]
%
% v = [ 1 1 1 ... 1 -1 -1 -1 ... -1 ]

% vectorization to fill gradient matrix
Gi = zeros(1,2*n);
Gj = zeros(1,2*n);
Gv = zeros(1,2*n);
for k = 1:n
    Gi(k) = k;
    Gi(k+n) = k;
    Gj(k) = k;
    Gj(k+n) = k+n;
    Gv(k) = 1;
    Gv(k+n) = -1;
end

% example: 5 x 10 = n x 2n
% gradient matrix
G = sparse(Gi,Gj,Gv);


% vectorization to fill 2n x 2n identity matrix
Ii_2n = zeros(1,2*n);
Ij_2n = zeros(1,2*n);
Iv_2n = zeros(1,2*n);
for k = 1:2*n
    Ii_2n(k) = k;
    Ij_2n(k) = k;
    Iv_2n(k) = 1;
end

% 2n x 2n identity matrix
I_2n = sparse(Ii_2n,Ij_2n,Iv_2n);


% matrix M = [G 0]^T * [G 0] + [I 0; 0 0]
M = sparse(Gi,Gj,Gv,n,2*n+n-1)'*sparse(Gi,Gj,Gv,n,2*n+n-1) + sparse(Ii_2n,Ij_2n,Iv_2n,2*n+n-1,2*n+n-1);

% discontinuity matrix
% i = 1:n-1
% D(i,i+1) = -1
% D(i,i+n) = 1
%
% i = [  1   2  3 ... n-1   1    2     3    ... n-1 ]
%
% j = [  2   3  4 ... n    n+1  n+2   n+3  ... 2n-1 ] 
%
% v = [ -1  -1 -1 ... -1    1     1     1   ...  1 ]

% vectorization to fill discontinuity matrix
Di = zeros(1,2*n-2);
Dj = zeros(1,2*n-2);
Dv = zeros(1,2*n-2);
for k = 1:n-1
    Di(k) = k;
    Di(k+n-1) = k;
    Dj(k) = k+1;
    Dj(k+n-1) = k+n;
    Dv(k) = -1;
    Dv(k+n-1) = 1;
end

% example: 4 x 10 = n-1 x 2n
% discontinuity matrix
D = sparse(Di,Dj,Dv,n-1,2*n);


% unkown variable vector    f       2n x 1
% discontinuity vector      Df      n-1 x 1

% b1 = -2 * g^T * I  <-- row ve    1 x m * m x n --> 1 x n
% b1' = [ b1 0]
% [ f Df  ] 

% b vector
temp_b = -2  * I_2n * g;
b = ones(1,2*n+n-1);
for i = 1:2*n
    b(i) = temp_b(i);
end 

% vectorization to fill n-1 x n-1 identity matrix
Ii_n_1 = zeros(1,n-1);
Ij_n_1 = zeros(1,n-1);
Iv_n_1 = zeros(1,n-1);
for k = 1:n-1
    Ii_n_1(k) = k;
    Ij_n_1(k) = k;
    Iv_n_1(k) = 1;
end

% n-1 x n-1 identity matrix
I_n_1 = sparse(Ii_n_1,Ij_n_1,Iv_n_1);

% x = quadprog(H,f,[],[],Aeq,beq)

f = quadprog(M,b',[-D -I_n_1; D -I_n_1],sparse(2*(n-1),1));

size(f)

% plot functions

    function v = format_vec(u)
        [m,~] = size(u);
        v = zeros(m/2,2);
        for i = 1:m/2
            v(i,1) = u(i);
            v(i,2) = u(i+m/2);
        end
    end 

f = f(1:2*n);

f = format_vec(f);

g = format_vec(g);

f1 = f(1:n,1);
f2 = f(1:n,2);

g1 = g(1:n,1);
g2 = g(1:n,2);

x_axis = 0:1/n:1;

% keep first n values
x_axis_temp = zeros(2*n,1);
for j = 1:n
    x_axis_temp(j) = x_axis(j);
end
% concatenate additional values
for j = n+1:2*n
    x_axis_temp(j) = x_axis(j-n+1);
end

x_axis = format_vec(x_axis_temp);

x1 = x_axis(1:n,1);
x2 = x_axis(1:n,2);

hold on

plot_g = plot([x1 x2], [g1 g2], 'Color', [0, 0.4470, 0.7410], 'DisplayName', 'g','LineWidth',1);

plot_f = plot([x1 x2], [f1 f2], 'Color', [0.8500, 0.3250, 0.0980], 'DisplayName', 'f','LineWidth',1);


title(['Denoising problem using quadprog results for ', num2str(n), ' sample points'])

hold off

end