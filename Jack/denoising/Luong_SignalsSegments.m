%Requires communication and optimization toolbox
opengl software
clear all

%Prompts user to input how many stepsizes should discretize [0,1]
stepsize = input("Enter the amount of steps the function will take: ");
points = stepsize + 1;
%Only doing square function for now because
signal = 0.5*(buildSquare(stepsize)+1);
%Adds white gaussian noise to the square function
signalWithError = addNoise(signal);
%Sets up and solves the l^2 optimization problem
alpha = 1;
beta = 1;
gamma = 100;
l2Solution = SolveL2Problem(signalWithError, stepsize, alpha, beta, gamma);
l2Signal = solutionRegularizer(l2Solution.recoveredSignal);
%l1Solution = SolveL1Problem(signalWithError, stepsize, alpha, beta, gamma);
%l1Signal = solutionRegularizer(l1Solution.recoveredSignal);

%Plot Line Segment Functions
hold on
axis equal
x = 0:1/stepsize:1;
for i = 1:stepsize
    if i == 1
        p1 = plot([x(i) x(i+1)], [signal(i, 1) signal(i,2)], 'Color', 'b', 'DisplayName', 'Original Signal', 'LineWidth', 1);
    else
        plot([x(i) x(i+1)], [signal(i, 1) signal(i,2)], 'Color', 'b', 'DisplayName', 'Original Signal', 'LineWidth', 1)
    end
end
%for i = 1:stepsize
%    if i == 1
%        p2 = plot([x(i) x(i+1)], [signalWithError(i, 1) signalWithError(i,2)], 'Color', 'r', 'DisplayName', 'Noisy Signal', 'LineWidth', 1);
%    else
%        plot([x(i) x(i+1)], [signalWithError(i, 1) signalWithError(i,2)], 'Color', 'r', 'DisplayName', 'Noisy Signal', 'LineWidth', 1)
%    end
%end
for i = 1:stepsize
    if i == 1
        p3 = plot([x(i) x(i+1)], [l2Signal(i, 1) l2Signal(i,2)], 'Color', 'g', 'DisplayName', 'l2 Signal', 'LineWidth', 1);
    else
        plot([x(i) x(i+1)], [l2Signal(i, 1) l2Signal(i,2)], 'Color', 'g', 'DisplayName', 'l2 Signal', 'LineWidth', 1)
    end
end
%for i = 1:stepsize
%    if i == 1
%        p4 = plot([x(i) x(i+1)], [l1Signal(i, 1) l1Signal(i,2)], 'Color', 'm', 'DisplayName', 'l1 Signal', 'LineWidth', 1);
%    else
%        plot([x(i) x(i+1)], [l1Signal(i, 1) l1Signal(i,2)], 'Color', 'm', 'DisplayName', 'l1 Signal', 'LineWidth', 1)
%    end
%end
hold off
title(['Denoising Algorithm with Line Segment Functions with parameters \alpha = ', num2str(alpha), ', \beta = ', num2str(beta), ', \gamma = ', num2str(gamma)])
legend([p1 p3])

%Plot Gradients of Denosied Functions
figure
hold on
scatter(1:stepsize, (computeDerivative(stepsize)*l2Solution.recoveredSignal)', 'LineWidth', 1)
%scatter(1:stepsize, (computeDerivative(stepsize)*l1Solution.recoveredSignal)', 'LineWidth', 1)
legend('Gradient vector of Recovered Signal with l2 term', 'Gradient vector of Recovered Signal with l1 term')
title('Gradient of Denoised Functions')
hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Functions

%Builds a square function using line segment functions
function squareFunction = buildSquare(points)
    squareFunction = zeros(points, 2);
    for i = 1:points
        if i < floor(points/2)
            squareFunction(i,1) = 1;
            squareFunction(i,2) = 1;
        elseif i == floor(points/2)
            squareFunction(i,1) = 1;
            squareFunction(i,2) = 1;
        else
            squareFunction(i,1) = -1;
            squareFunction(i,2) = -1;
        end
    end
end

%Builds an n x 2n derivative matrix for [startpoint endpoint] functions
function gradientMatrix = computeDerivative(stepsize)
    gradientMatrix = zeros(stepsize, 2*stepsize);
    for i = 1:stepsize
        gradientMatrix(i, i) = -1;
        gradientMatrix(i, stepsize + i) = 1;
    end
end

%Builds a n x 2n matrix to evaluate endpoint differences
function endpointMatrix = createEndpointMatrix(stepsize)
    endpointMatrix = zeros(stepsize, 2*stepsize);
    for i = 2:stepsize
        endpointMatrix(i,i) = -1;
        endpointMatrix(i, stepsize + i - 1) = 1;
    end
end

%Adds Gaussian noise to a function, preserves continuity of orignal
%function
function noisyFunction = addNoise(originalFunction)
    [m,n] = size(originalFunction);
    noisyFunction = zeros(m,n);
    for i = 1:m
        for j = 1:n
            if i == 1 && j == 1
                noisyFunction(i,j) = originalFunction(i,j) + normrnd(0,0.1);
            elseif i == m && j == n
                noisyFunction(i,j) = originalFunction(i,j) + normrnd(0,0.1);
            elseif j == n
                noise = normrnd(0,0.1);
                noisyFunction(i,j) = originalFunction(i,j) + noise;
                noisyFunction(i+1,j-1) = originalFunction(i+1,j-1) + noise;
            else    
            end
        end
    end
end

%Sets up, and solves the l2 optimization problem.  Returns the optimization
%solution
function l2Solution = SolveL2Problem(noisyFunction, stepsize, weight1, weight2, weight3)
    noisyFunctionReg = zeros(2*stepsize, 1);
    for i = 1:stepsize
        noisyFunctionReg(i) = noisyFunction(i,1);
        noisyFunctionReg(stepsize + i) = noisyFunction(i,2); 
    end
    l2problem = optimproblem('ObjectiveSense', 'min');
    recoveredSignal = optimvar('recoveredSignal', 2*stepsize, 1);
    endpointPenalizer = optimvar('errorPenalizer', stepsize, 1);
    %Set up terms used in objective
    error1 = recoveredSignal - noisyFunctionReg;
    error2 = computeDerivative(stepsize)*recoveredSignal;
    error3 = createEndpointMatrix(stepsize)*recoveredSignal;
    %Constraints
    endpointConstraint = optimconstr(stepsize - 1);
    endpointPenalty = optimconstr(2*stepsize);
    for i = 1:stepsize
        if i > 1
            endpointConstraint(i) = recoveredSignal(i) == recoveredSignal(stepsize + i - 1);
        end
    end
    for i = 1:2:2*stepsize
        endpointPenalty(i) = error3((i+1)/2) <= endpointPenalizer((i+1)/2);
        endpointPenalty(i+1) = -error3((i+1)/2) <= endpointPenalizer((i+1)/2);
    end
    l2problem.Objective = weight1*sum(error1.^2) + weight2*sum(error2.^2) + weight3*sum(endpointPenalizer);
    %l2problem.Constraints.cons1 = endpointConstraint; %Toggles continuity
    l2problem.Constraints.cons2 = endpointPenalty; %Encourages continuity
    l2Solution = solve(l2problem);
end

%Sets up, and solves the l1 optimization problem.  Returns the optimization
%solution
function l1Solution = SolveL1Problem(noisyFunction, stepsize, weight1, weight2, weight3)
    noisyFunctionReg = zeros(2*stepsize, 1);
    for i = 1:stepsize
        noisyFunctionReg(i) = noisyFunction(i,1);
        noisyFunctionReg(stepsize + i) = noisyFunction(i,2); 
    end
    l1problem = optimproblem('ObjectiveSense', 'min');
    recoveredSignal = optimvar('recoveredSignal', 2*stepsize, 1);
    recoveredSignalDeriv = optimvar('recoveredSignalDeriv', stepsize, 1);
    endpointPenalizer = optimvar('errorPenalizer', stepsize, 1);
    %Set up terms used in objective
    error1 = recoveredSignal - noisyFunctionReg;
    error2 = computeDerivative(stepsize)*recoveredSignal;
    error3 = createEndpointMatrix(stepsize)*recoveredSignal;
    %Constraints
    endpointConstraint = optimconstr(stepsize - 1);
    absoluteValueConstraint = optimconstr(2*stepsize);
    endpointPenalty = optimconstr(2*stepsize);
    for i = 1:stepsize
        if i > 1
            endpointConstraint(i) = recoveredSignal(i) == recoveredSignal(stepsize + i - 1);
        end
    end
    for i = 1:2:2*stepsize
        absoluteValueConstraint(i) = -error2((i+1)/2) <= recoveredSignalDeriv((i+1)/2);
        absoluteValueConstraint(i+1) = error2((i+1)/2) <= recoveredSignalDeriv((i+1)/2);
    end
    for i = 1:2:2*stepsize
        endpointPenalty(i) = error3((i+1)/2) <= endpointPenalizer((i+1)/2);
        endpointPenalty(i+1) = -error3((i+1)/2) <= endpointPenalizer((i+1)/2);
    end
    l1problem.Objective = weight1*sum(error1.^2) + weight2*sum(recoveredSignalDeriv) + weight3*sum(endpointPenalizer);
    %l1problem.Constraints.cons1 = endpointConstraint; %Toggles continuity
    l1problem.Constraints.cons2 = absoluteValueConstraint; %For l1 derivative
    l1problem.Constraints.cons3 = endpointPenalty; %Encourages continuity
    l1Solution = solve(l1problem);
end

%Converts [startpoint endpoint] 2n vector to n x 2 [startpoind; endpoint]
%vector
function convertedVector = solutionRegularizer(vector)
    [n, ~] = size(vector);
    convertedVector = zeros(n/2, 2);
    for i = 1:n/2
        convertedVector(i,1) = vector(i);
        convertedVector(i,2) = vector(i + n/2);
    end
end