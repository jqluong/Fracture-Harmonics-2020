function f_value = prob_noise_min_seg(stepsize)

% vectorization of segment vertices: (assuming endpoints have a single vetex, mid points two vertices)
%
% [ f_12 f_21 f_23 f_32 ... f_ij f_ji f_jk f_kj f_kl f_lk ... f_nm f_mn ]
%
% possible constraint f_ij = f_ji for all i = 2:n-1, j = 1:m
%

% number of vertex (regular)
n0 = (1/stepsize) + 1;

% number of vertex (segmented)
n = 2*(n0 - 1);

% create noisy input 
% first use standard square with desired amplitude and period
x = 0:2*pi/(n0-1):2*pi;
y = 5 + square(x);
% then add random noise
g0 = y + (0.1) * randn(1,n0);
% create multivalued-segmented function using noisy input
g = zeros(1,n);
g(1) = g0(1);
for i = 2:n0-1
    g(2*i-2) = g0(i);
    g(2*i-1) = g(2*i-2);
end
g(n) = g0(n0);
        
% create matrix to take the discrete L1 of "continuity between segments"
discL1 = zeros(n0,n);
    for j = 2:2:n0-1
        discL1(j,j) = 1;
        discL1(j,j+1) = -1;
    end
    
% create optimization problems:
% for L2 term energy expression
minL2prob = optimproblem('ObjectiveSense','minimize');
% for L1 term energy expression
minL1prob  = optimproblem('ObjectiveSense','minimize');

% create optimization variables
f = optimvar('f',n,1);                  
Y = optimvar('Y',n0,1); 

J1 = f-transpose(g);                    

J2 = discL1*f;

J2star = discL1*f;

% define objective functions:
% for L2 term problem
minL2prob.Objective = sum(J1.^2) + sum(J2.^2);
% for L1 term problem
minL1prob.Objective = sum(J1.^2) + sum(Y);

% constraint for absolute value in L1
absConstr = optimconstr(2*n0);
for i = 1:2:2*n0-1
    j = (i+1)/2;
    absConstr(i) = -J2star(j)  <= Y(j);
    absConstr(i+1) = J2star(j) <= Y(j);
end
minL1prob.Constraints.absL1Constr = absConstr;

% constraint multiple edpoints at same sample pt coincide
sameConstr = optimconstr(n0);
for j = 1:n0
    sameConstr(j) = J2star(j) == 0;
end

% solve min problems:
solL2min = solve(minL2prob);
solL1min = solve(minL1prob);

% function used to format vector 
function formated_vect = format_vector(u)
[p,m] = size(u);
alpha_formated = zeros(p/2,2);
for k = 1:p/2
    alpha_formated(k,1) = u(2*k-1);
    alpha_formated(k,2) = u(2*k);
end
formated_vect = alpha_formated;
end
    
% format L2 term and L1 term vectors
fin_g = format_vector(transpose(g));
finSolL2 = format_vector(solL2min.f);
finSolL1 = format_vector(solL1min.f);

%check result (use for debug)
%fin_g
%finSolL2
%finSolL1

% plot results
x_axis = 0:1/(n-1):1;
hold on
% noisy funtion
for i = 1:n/2
    if i == 1
        plot_noisy = plot([x_axis(i) x_axis(i+1)], [fin_g(i, 1) fin_g(i,2)], 'Color', [0, 0.4470, 0.7410], 'DisplayName', 'noisy function');
    else
        plot([x_axis(i) x_axis(i+1)], [fin_g(i, 1) fin_g(i,2)], 'Color', [0, 0.4470, 0.7410], 'DisplayName', 'noisy function');
    end
end
% L2 term min function
for i = 1:n/2
    if i == 1
        plot_L2 = plot([x_axis(i) x_axis(i+1)], [finSolL2(i, 1) finSolL2(i,2)], 'Color', [0.8500, 0.3250, 0.0980], 'DisplayName', '(f prime)^2 unknown function');
    else
        plot([x_axis(i) x_axis(i+1)], [finSolL2(i, 1) finSolL2(i,2)], 'Color', [0.8500, 0.3250, 0.0980], 'DisplayName', '(f prime)^2 unknown function');
    end
end
% L1 term min function
for i = 1:n/2
    if i == 1
        plot_L1 = plot([x_axis(i) x_axis(i+1)], [finSolL1(i, 1) finSolL1(i,2)], 'Color', [0.8500, 0.3250, 0.0980], 'DisplayName', '|f prime|^2 unknown function');
    else
        plot([x_axis(i) x_axis(i+1)], [finSolL1(i, 1) finSolL1(i,2)], 'Color', [0.8500, 0.3250, 0.0980], 'DisplayName', '|f prime| unknown function');
    end
end
hold off       
legend([plot_noisy plot_L2 plot_L1])        

end