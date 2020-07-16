%
% function laplacian = discrete_laplacian(f,x)
%
% This function takes the discrete laplacian of a function f using the
% finite difference method
%
% paramaters: f - array of values f_n i.e. values of f at each x_n
%             x - array of values x_n that determine step size
%
% return: laplacian - column vector with values of discrete laplacian of f
%
% example: 
%
%   >> x = [0,0.5,1,1.5,2,2.5];
%   >> y = zeros(length(x),1);
%   >> for i = 1:length(x)
%          y(i,1)=(x(1,i))^2;
%      end
%   >> discrete_laplacian(y,x);
%   
%   laplacian =
%       -1
%        2
%        2
%        2
%        2
%       -9
%         

function laplacian = discrete_laplacian(f,x)

n = length(f);      % determine lenght of vector input
L = zeros(n);       % create matrix L 
z = zeros(n-1,1);   % create vector for "middle points" step

% fill "middle points" step vector
for k = 2:n-1
    z(k) = ( (x(k) + (x(k+1) - x(k))/2) ) - ( x(k) - ((x(k) - x(k-1))/2) );
end
for k = 1
    z(k) = ((x(k+1) - x(k))/2) + ((x(k+1) - x(k))/2) ;
end
for k = n
    z(k) =  ((x(k) - x(k-1))/2) + ((x(k) - x(k-1))/2);
end

% fill matrix L
for i = 2:n-1
    L(i,i-1) = 1/( ( x(i) - x(i-1) )*z(i) ); 
    L(i,i+1) = 1/( ( x(i+1) - x(i) )*z(i) ); 
    L(i,i) = - L(i,i+1)  - L(i,i-1);
end
for i = 1
    L(i,i)= 1/( ( x(i+1) - x(i) )*z(i) );
    L(i,i+1) = -1/( ( x(i+1) - x(i) )*z(i) );
end
for i=n
    L(i,i-1) = 1/( ( x(i) - x(i-1) )*z(i) );
    L(i,i) = -1/( ( x(i) - x(i-1) )*z(i) );
end

% return discrete Laplacian
laplacian = L*f

end

