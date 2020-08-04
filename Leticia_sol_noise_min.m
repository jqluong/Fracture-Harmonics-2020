function f_values = sol_noise_min (stepsize)

% number of vertex
n = (1/stepsize) + 1;

% define optimization options
options = optimoptions(@fminunc,'GradObj','off','MaxIter',1e+4);       % GradObj defines gradient of obj funtion

% objective function used in min
    function [funL2] = objFunL2(u,g)                        % note:  u is dummy variable

        funL2 = (norm(u-g))^2 + (norm(gradient(u)))^2;      % objective: ||f-g||^2 + (f')^2
        
      % grad = 2*(u - g) + ( 2*del2(u) );                   % gradient:  f - g + (1/2)f"
        
    end

    function [funL1] = objFunL1(v,g)
        
        funL1 = (norm(v-g))^2 + sum(abs(gradient(v)));      % objective: ||f-g||^2 + |f'|^2
        
        % let's not try gradient...
        
    end

% create noisy function input:
%
% first use standard square with desired amplitude and period
x = 0:2*pi/(n-1):2*pi;
y = 5 + (square(x));
%
% then add random noise
g = y + (0.1) * randn(1,n);
%

% create function handle to pass extra arguments
fL2 = @(u)objFunL2(u,g);

fL1 = @(v)objFunL1(v,g);

% initialize starting function
f0 = g;

% find min
fL2Min = fminunc(fL2,f0,options);

fL1Min = fminunc(fL1,f0,options);

x_axis =  0:1/(n-1):1;
plot (x_axis,g);
hold on
plot (x_axis,fL2Min);
plot (x_axis,fL1Min);
legend('noisy function','(f prime)^2 unknown function', '|f prime| unknown function');
hold off
   
end