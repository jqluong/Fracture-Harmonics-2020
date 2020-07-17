a = 0;  %left endpoint
b = 2*pi;  %right endpoint
n = 100; %number of partitions
h = (b-a)/n;
x = linspace(a,b,n+1)'; 
fx = zeros(1,length(x))';   %f(x) = sin(x)


for i = 1:length(x)
    fx(i) = f(x(i));
end

fx_without_endpoints = fx(2:n);  %exclude endpoints at which the second derivative is not defined

D = (1/(h^2))*(diag(-2*ones(1,n-1)) + diag(ones(1,n-2),1) + diag(ones(1,n-2),-1)); %matrix to compute second derivatives using centered difference approx

second_derivative = D*fx_without_endpoints; %Matrix multiplication that computes a second derivative vector


%% Plot graph of f(x) = sin(x) and f''(x) = -sin(x)

plot(x(2:n),fx_without_endpoints);

hold on

plot(x(2:n),second_derivative);




 function y = f(x)
     y = sin(x);
 end