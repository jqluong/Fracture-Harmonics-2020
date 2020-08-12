function f_value = penalize_seg_discrep(stepsize, a, b, c)

% vectorization of segment vertices: (assuming endpoints have a single vetex, mid points two vertices)
%
% [ f_12 f_21 f_23 f_32 ... f_ij f_ji f_jk f_kj f_kl f_lk ... f_nm f_mn ]
%
% possible constraint f_ij = f_ji for all i = 2:n-1, j = 1:m
%

% new parameters
tau = a;
lambda = b;
eta = c;

% number of vertex (regular)
n0 = (1/stepsize) + 1;

% number of vertex (segmented)
n = 2*(n0 - 1);

% create noisy input 
% first use standard square with desired amplitude and period
x = 0:2*pi/(n0-1):2*pi;
y = 0.5+ 0.25*square(x);
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
discL1 = zeros(n0-1,n);
    for j = 1:n0-1
        discL1(j,2*j-1) = 1;
        discL1(j,2*j) = -1;
    end


% crate matrix to penalize discrepancy
S = zeros(n0,n);
for i = 2:n0-1
    S(i,2*i-2) = 1 ;
    S(i,2*i-1) = -1;
end 

    
% create optimization problems:
% for L2 term energy expression
minL2prob = optimproblem('ObjectiveSense','minimize');
% for L1 term energy expression
minL1prob  = optimproblem('ObjectiveSense','minimize');

% create optimization variables
f = optimvar('f',n,1);                  
Y = optimvar('Y',n0-1,1); 
W = optimvar('W',n0-2,1); 

J1 = f-transpose(g);  

J2 = discL1*f;  % vector (n0-1) length

J2star = discL1*f;

J3 = S*f;

% define objective functions:
% for L2 term problem
minL2prob.Objective = tau*sum(J1.^2) + lambda*sum(J2.^2) + eta*sum(W);
% for L1 term problem
minL1prob.Objective = tau*sum(J1.^2) + lambda*sum(Y) + eta*sum(W);

% constraint for absolute value in L1
absConstr = optimconstr(2*(n0-1));
for i = 1:2:2*n0-3                          %2*(n0-1)
    j = (i+1)/2;                            % (2*n0 - 3 +1)/2 = n0-1  % 2*n0 - 3 = 2*(n0-1) -1
    absConstr(i) = -J2star(j)  <= Y(j);
    absConstr(i+1) = J2star(j) <= Y(j);     % j = n0-1 on last iter
end
minL1prob.Constraints.absL1Constr = absConstr;

% constraint for absolute value in L1
penaltyConstr = optimconstr(2*(n0-2));
for i = 1:2:2*n0-5                      
    j = (i+1)/2;                            % (2*n0 - 4 -1) = 2*n0 -5
    penaltyConstr(i) = -J3(j)  <= W(j);
    penaltyConstr(i+1) = J3(j) <= W(j);     % j = n0-1 on last iter
end
minL1prob.Constraints.penaltyL1Constr = penaltyConstr;
minL2prob.Constraints.penaltyL2Constr = penaltyConstr;

% constraint multiple edpoints at same sample pt coincide
sameConstr = optimconstr(n/2-1);  
for j = 1:n/2-1
   %sameConstr(j) = J2star(j) == 0;     % replaced to avoid issues
   sameConstr(j) = f(2*j) == f(2*j+1);
   
   
   % sameConstr(1) = f(2)   == f(3)
   % sameConstr(2) = f(4)   == f(5)
   % sameConstr(3) = f(6)   == f(7)
   % sameConstr(4) = f(8)   == f(9)
   % sameConstr(i) = f(2*i) == f(2*i+1)
   % sameConstr(n/2-1) = f(n-2) == f(n-1)
   
   % 2*i+1 = n-1 -> 2*i = n-2 -> i = n/2 -1
   
end

% remove constraint if desired:
%minL2prob.Constraints.sameConstr = sameConstr;
%minL1prob.Constraints.sameConstr = sameConstr;


% solve min problems:
solL2min = solve(minL2prob);
solL1min = solve(minL1prob);

% function used to format vector 
function formated_vect = format_vector(u)
[p,~] = size(u);
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


% plot results
x_axis = 0:2/(n-1):2;
hold on
% noisy funtion
for i = 1:n/2
    if i == 1
        plot_noisy = plot([x_axis(i) x_axis(i+1)], [fin_g(i, 1) fin_g(i,2)], 'Color', [0, 0.4470, 0.7410], 'DisplayName', 'noisy function','LineWidth',1.5);
    else
        plot([x_axis(i) x_axis(i+1)], [fin_g(i, 1) fin_g(i,2)], 'Color', [0, 0.4470, 0.7410], 'DisplayName', 'noisy function','LineWidth',1.5);
    end
end
title(['Denoising problem results for ', num2str(n0-1) ,' sample points with parameters \tau = ', num2str(tau), ', \lambda = ', num2str(lambda), ', \eta = ', num2str(eta)])
legend(plot_noisy)
figure;
hold off;
% L2 term min function
hold on;
for i = 1:n/2
    if i == 1
        plot_L2 = plot([x_axis(i) x_axis(i+1)], [finSolL2(i, 1) finSolL2(i,2)], 'Color', [0.8500, 0.3250, 0.0980], 'DisplayName', '(f prime)^2 unknown function','LineWidth',1.5);
    else
        plot([x_axis(i) x_axis(i+1)], [finSolL2(i, 1) finSolL2(i,2)], 'Color', [0.8500, 0.3250, 0.0980], 'DisplayName', '(f prime)^2 unknown function','LineWidth',1.5);
    end
end
title(['Denoising problem results for ', num2str(n0-1) ,' sample points with parameters \tau = ', num2str(tau), ', \lambda = ', num2str(lambda), ', \eta = ', num2str(eta)])
legend(plot_L2)
figure;
hold off;
% L1 term min function
hold on;
for i = 1:n/2
    if i == 1
        plot_L1 = plot([x_axis(i) x_axis(i+1)], [finSolL1(i, 1) finSolL1(i,2)], 'Color', [0.4940, 0.1840, 0.5560], 'DisplayName', '|f prime|^2 unknown function','LineWidth',1.5);
    else
        plot([x_axis(i) x_axis(i+1)], [finSolL1(i, 1) finSolL1(i,2)], 'Color', [0.4940, 0.1840, 0.5560], 'DisplayName', '|f prime| unknown function','LineWidth',1.5);
    end
end
title(['Denoising problem results for ', num2str(n0-1) ,' sample points with parameters \tau = ', num2str(tau), ', \lambda = ', num2str(lambda), ', \eta = ', num2str(eta)])
legend(plot_L1)
hold off       


end