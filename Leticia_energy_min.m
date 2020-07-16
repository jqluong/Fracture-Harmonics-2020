%
% function f_values = energy_min(a,b,x)
%
% This function solves the laplace equation Δf=0 using energy minizer
% method
%
% paramaters: a - boundary condition at x_0
%             b - boundary condition at x_n
%             x - array of values x_n that determine step size
%
% return: f_values - column vector with values of function f at each x_n
%
% example: 
%
%   >> x = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
%   >> energy_min(0,1,x);
%
%   f_values =
%
%            0
%          0.1
%          0.2
%          0.3
%          0.4
%          0.5
%          0.6
%          0.7
%          0.8
%          0.9
%            1
%

function f_values = energy_min (a,b,x)

n = length(x);          % determine lenght of vector input
G = zeros(n-1,n);       % create matrix G
A = zeros(n-1);         % create matrix A    
f_0 = zeros(n,1);       % create vector [f1 0 ... 0 fn]
f = zeros(n,1);         % create vector [f1 f2 ... fn-1 fn]


% boundary conditions
f_0(1) = a;
f_0(n) = b;

% fill matrix G
for i = 1:n-1
    G(i,i) = 1/(x(i+1) - x(i));
    G(i,i+1) = -1/(x(i+1) - x(i));
end

% fill matrix A
for i = 1:n-1
    A(i,i) = x(i+1) - x(i);
end

Q = transpose(G)*A*G;   % define matrix Q to simplify

z = Q*f_0;              % define vector z to simplify

% (1/2)*f_unk ^T * Q(2:n-1,2:n-1) * f_unk + z(2:n-1)*f_unk
%
% Q(2:n-1,2:n-1) * f_unk + z(2:n-1) =(set)= 0

f_unk = - Q(2:n-1,2:n-1) \ z(2:n-1,1);  % solve linear system

% fill vector with function values
f(1) = f_0(1); 
f(n) = f_0(n);
for i = 2:n-1
    f(i)=f_unk(i-1);
end

%  return vector with function values
f_values = f

end