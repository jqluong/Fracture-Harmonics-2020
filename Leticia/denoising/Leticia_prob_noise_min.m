function f_value = prob_noise_min(stepsize)

% number of vertex
n = (1/stepsize) + 1;

% create gradient funtion that takes optimvar
G = zeros(n,n);
for j = 2:n-1
   G(j,j+1) = 0.5;
   G(j,j) = 0;
   G(j,j-1) = -0.5;
end
   G(1,1) = -1;
   G(1,2) = 1;
   G(n,n) = 1;
   G(n,n-1) = -1;

% create noisy input 
% first use standard square with desired amplitude and period
x = 0:2*pi/(n-1):2*pi;
y = 5 + square(x);
% then add random noise
g = y + (0.1) * randn(1,n);

% create optimization problems:
% for L2 term energy expression
minL2prob = optimproblem('ObjectiveSense','minimize');
% for L1 term energy expression
minL1prob  = optimproblem('ObjectiveSense','minimize');

% create optimization variables
f = optimvar('f',n,1);                  
Y = optimvar('Y',n,1); 

J1 = f-transpose(g);                    

J2 = G*f;

J2star = G*f;

% define objective functions:
% for L2 term problem
minL2prob.Objective = sum(J1.^2) + sum(J2.^2);
% for L1 term problem
minL1prob.Objective = sum(J1.^2) + sum(Y);

% constraint for absolute value in L1
absConstr = optimconstr(2*n);
for i = 1:2:2*n-1
    j = (i+1)/2;
    absConstr(i) = -J2star(j)  <= Y(j);
    absConstr(i+1) = J2star(j) <= Y(j);
end
minL1prob.Constraints.absL1Constr = absConstr;

% solve min problems:
solL2min = solve(minL2prob);
solL1min = solve(minL1prob);

% plot results
x_axis =  0:1/(n-1):1;
plot (x_axis,g);
hold on
plot (x_axis,solL2min.f);
plot (x_axis,solL1min.f);
legend('noisy function','(f prime)^2 unknown function', '|f prime| unknown function');
hold off

end