function f_values = noise_grad_desc(stepsize)

% suggested values for parameters
 tau_1 = 0.1;
 lamb_1 = 1.8;
 tau_2 = 0.1;
 lamb_2 = 0.1;
 e = 0.001;
 iterations = 100;

% number of nodes
n = (1/stepsize) + 1;

% used for l1 norm
grad_H2 = zeros(1,n);

% create noisy input 
% first use standard square with desired amplitude and period
x = 0:2*pi/(n-1):2*pi;
y = 5 + (square(x));
% then add random noise
g = y + (0.1) * randn(1,n);


% case min ||f-g||^2 + (f')^2 

% initialize unknown function to noisy function
f = g;

% iterate gradient descent
for i = 1:iterations
    
    f_old = f;
    
    % update f
    f = f_old - tau_1*((f_old-g)-lamb_1*4*del2(f_old));  
    
    J = (norm(f_old-g))^2 + (norm(gradient(f_old)))^2;
end

% case min ||f-g||^2 + |f'| 

% initialize unknown function to noisy function
h = g;

% iterate gradient descent
for j = 1:iterations
    
    h_old = h;
    
    % grad h vector
    grad_h = gradient(h_old);
    
    % u vector to be used in grad_H2
    u = zeros(n,1);
    for m = 1:n
        u(m) = grad_h(m)/sqrt(e^2+(grad_h(m))^2);
    end
    
    
    % grad_H2 vector for gradient descent
    grad_H2 = gradient(u);
    
    
    % update h
    h = h_old - tau_2*( (h_old-g) - lamb_2*( transpose(grad_H2) ) ) ;
    
    H = (norm(h_old-g))^2 + sum(abs(grad_h));
    
end

% plot functions
x_axis =  0:1/(n-1):1;
plot (x_axis,g);
hold on
plot (x_axis,f);
plot (x_axis,h);
legend('noisy function','(f prime)^2 unknown function','|f prime| unknown function');
hold off
    
end