%
% function f_values = laplace_eq(a,b,x)
%
% This function solves the laplace equation Δf=0 using finite difference
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
%   >> laplace_eq(0,1,x);
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

function f_values = laplace_eq(a,b,x)

n = length(x);              % determine lenght of vector input
y = zeros(n-2,1);           % will solve M*y = 0
f = zeros(n-1,1);           % vector [f1 f2 ... fn]
M = zeros(n-2);             % create matrix M
L = zeros(n);               % create matrix L 

% boundary conditions       
f(1) = a;
f(n) = b;

% adjust y to boundary conditions
y(1)=-a;    
y(n-2)=-b;

% fill matrix L
for i = 2:n-1
    L(i,i-1) = ( lcm(sym(x(i) - x(i-1)),sym(x(i+1) - x(i))) ) / ( x(i) - x(i-1) ) ; 
    L(i,i+1) = ( lcm(sym(x(i) - x(i-1)),sym(x(i+1) - x(i))) ) / ( x(i+1) - x(i) ) ; 
    L(i,i) = - L(i,i+1)  - L(i,i-1);
end

% fill matrix M
for i = 1:n-2
    for j = 1:n-2
    M(i,j) = L(i+1,j+1);
    end
end

% solve linear system, z is dummy variable
z = M \ y;

% fill vector [f1 f2 ... fn] with f2 ... fn-1
for i = 2:n-1
    f(i) = z(i-1);
end

% return vector [f1 f2 ... fn]
f_values = f

end
