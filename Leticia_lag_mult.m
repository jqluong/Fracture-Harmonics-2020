%
% function f_values = lag_mult(x)
%
% This function solves the laplace equation Δf=0 subject to ||f||=1 using lagrange multipliers
%
% paramaters: x - array of values x_n that determine step size
%
% return: f_values - column vector with values of function f at each x_n
%
% example: 
%
%   >> x = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
%   >> lag_mult(0,1,x);
%
%   f_values =
%
%             1
%             1
%             1
%             1
%             1
%             1
%             1
%             1
%             1
%             1
%             1
%

function f_values = lag_mult(x)

n = length(x);          % determine lenght of vector input
G = zeros(n-1,n);       % create matrix G
A = zeros(n-1);         % create matrix A   
f = zeros(n,1);         % create vector [f1 f2 ... fn-1 fn]
h = (x(n)-x(1))/n       % step size if assuming equal step sizes


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

% minimize E(f) = 0.5*f^T*Q*f
%
% subject to g(f) = f^T*f-1 = 0
%
% Lamb = 0.5*f^T*Q*f - lamb^T(f^T*f-1);    % Lagrange function
%
%Q*f-lamb^T*f = 0;
%
%f^T*f-1 = 0;

% find eigenvalues and eigenvectors of matrix Q
[V,D] = eig(Q)

% create vector of eigenvalues
lamb = diag(D);

% determine and locate the minimum of eigenvalues
[lamb_min,index] = min(lamb);

% return the corresponding eigenvector of minimum eigenvalue 
f = (1/sqrt(h))*V(:,index)

% plot first five eigenvectors
for i = 1:5
    plot (x,(1/sqrt(h))*V(:,i));
    hold on
end

legend('eigen f_1','eigen f_2','eigen f_3','eigen f_4','eigen f_5')

end
