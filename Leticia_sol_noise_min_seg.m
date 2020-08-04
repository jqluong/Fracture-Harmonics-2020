function f_values = sol_noise_min_seg(stepsize)

% number of vertex
n = (1/stepsize) + 1;

% [ f_12 f_21 f_23 f_32 ... f_ij f_ji f_jk f_kj f_kl f_lk ... f_nm f_mn ]
%
% possible constraint f_ij = f_ji for all i = 2:n-1, j = 1:m
%

% create matrix to take the discrete L1 of "continuity between segments"
G = zeros(n,2*n);
for k = 2:n-1
G(k,k) = 1;
G(k,k+1) = -1;
end

% create noisy function input:
%
% first use standard square with desired amplitude and period
x = 0:2*pi/(n-1):2*pi;
y = 5 + (square(x));
%
% then add random noise
g0 = y + (0.1) * randn(1,n);
%
% create multivalued-segmented function using noisy input
g = zeros(1,2*n);
g(1) = g0(1);
for i = 2:n-1
    g(2*i-2) = g0(i);
    g(2*i-1) = g(2*i-2);
end
        g(2*n) = g0(n);
              
% initialize f
  f0  = transpose(g);
        
% objective function used in min
function [funL2] = objFunL2(u,g)
    [p,m] = size(g);
    
    q = m/2;
         
% create matrix to take the discrete L1 of "continuity between segments"
Grad = zeros(q,m);
    for j = 2:q-1
        Grad(j,j) = 1;
        Grad(j,j+1) = -1;
    end
        
   
grad_u = Grad*u;
        
funL2 = (norm(u-g))^2 + (norm(grad_u))^2 ;      % objective: ||f-g||^2 + (f')^2
        
end
                
% constraints
A = [];
b = [];
Aeq = G;
beq = zeros(n,1);
%
        
        
% create function handle to pass extra arguments
fL2 = @(u)objFunL2(u,g);

%Obj0 = @(c)0;

%fL1 = @(v)objFunL1(v,g);   % suggested fix of "Converged to an infeasible point"
        
%f0 = fmincon(Obj0,f0,A,b,Aeq,beq);
        
% minimization
fL2Min = fmincon(fL2,transpose(g),A,b,Aeq,beq);
        
        
        
% function used to format vector 
function formated_vect = format_vector(alpha)
[w,z] = size(alpha);
r = w/2;
alpha_formated = zeros(r,2);
for j = 1:r
    alpha_formated(j,1) = alpha(2*j-1);
    alpha_formated(j,2) = alpha(2*j);
end
formated_vect = alpha_formated;
end
    
% format L2 term and L1 term vectors
L2fun = format_vector(fL2Min);
        
gfun = format_vector(transpose(g));
        
        
x_axis = 0:1/(n-1):1;
hold on
% noisy funtion
for i = 1:n-1
    if i == 1
        plot_noisy = plot([x_axis(i) x_axis(i+1)], [gfun(i, 1) gfun(i,2)], 'Color', [0, 0.4470, 0.7410], 'DisplayName', 'noisy function')
    else
        plot([x_axis(i) x_axis(i+1)], [gfun(i, 1) gfun(i,2)], 'Color', [0, 0.4470, 0.7410], 'DisplayName', 'noisy function')
    end
end
% L2 term min function
for i = 1:n-1
    if i == 1
        plot_L2 = plot([x_axis(i) x_axis(i+1)], [L2fun(i, 1) L2fun(i,2)], 'Color', [0.8500, 0.3250, 0.0980], 'DisplayName', '(f prime)^2 unknown function')
    else
        plot([x_axis(i) x_axis(i+1)], [L2fun(i, 1) L2fun(i,2)], 'Color', [0.8500, 0.3250, 0.0980], 'DisplayName', '(f prime)^2 unknown function')
    end
end
% L1 term min function
hold off       
legend([plot_noisy plot_L2])        

end