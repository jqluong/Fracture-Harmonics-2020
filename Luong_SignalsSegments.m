%Requires communication and optimization toolbox
opengl software
clear all

%Prompts user to input how many stepsizes should discretize [0,1]
stepsize = input("Enter the amount of steps the function will take: ");
points = stepsize + 1;
%Only doing square function for now because
signal = buildSquare(stepsize);
%Adds white gaussian noise to the square function
signalWithError = addNoise(signal);
%Sets up and solves the l^2 optimization problem
l2Solution = SolveL2Problem(signalWithError, stepsize);
l2Signal = solutionRegularizer(l2Solution.recoveredSignal);
l1Solution = SolveL1Problem(signalWithError, stepsize);
l1Signal = solutionRegularizer(l1Solution.recoveredSignal);

%Plot Line Segment Functions
hold on
x = 0:1/stepsize:1;
for i = 1:stepsize
    if i == 1
        p1 = plot([x(i) x(i+1)], [signal(i, 1) signal(i,2)], 'Color', 'b', 'DisplayName', 'Original Signal')
    else
        plot([x(i) x(i+1)], [signal(i, 1) signal(i,2)], 'Color', 'b', 'DisplayName', 'Original Signal')
    end
end
%for i = 1:stepsize
%    if i == 1
%        p2 = plot([x(i) x(i+1)], [signalWithError(i, 1) signalWithError(i,2)], 'Color', 'r', 'DisplayName', 'Noisy Signal')
%    else
%        plot([x(i) x(i+1)], [signalWithError(i, 1) signalWithError(i,2)], 'Color', 'r', 'DisplayName', 'Noisy Signal')
%    end
%end
for i = 1:stepsize
    if i == 1
        p3 = plot([x(i) x(i+1)], [l2Signal(i, 1) l2Signal(i,2)], 'Color', 'g', 'DisplayName', 'l2 Signal')
    else
        plot([x(i) x(i+1)], [l2Signal(i, 1) l2Signal(i,2)], 'Color', 'g', 'DisplayName', 'l2 Signal')
    end
end
for i = 1:stepsize
    if i == 1
        p4 = plot([x(i) x(i+1)], [l1Signal(i, 1) l1Signal(i,2)], 'Color', 'm', 'DisplayName', 'l1 Signal')
    else
        plot([x(i) x(i+1)], [l1Signal(i, 1) l1Signal(i,2)], 'Color', 'm', 'DisplayName', 'l1 Signal')
    end
end
hold off
legend([p1 p3 p4])

%Functions
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

function gradientMatrix = computeDerivative(stepsize)
    gradientMatrix = zeros(stepsize, 2*stepsize);
    for i = 1:stepsize
        gradientMatrix(i, i) = -1;
        gradientMatrix(i, stepsize + i) = 1;
    end
end

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

function l2Solution = SolveL2Problem(noisyFunction, stepsize)
    noisyFunctionReg = zeros(2*stepsize, 1);
    for i = 1:stepsize
        noisyFunctionReg(i) = noisyFunction(i,1);
        noisyFunctionReg(stepsize + i) = noisyFunction(i,2); 
    end
    l2problem = optimproblem('ObjectiveSense', 'min');
    recoveredSignal = optimvar('recoveredSignal', 2*stepsize, 1);
    %Set up terms used in objective
    error1 = recoveredSignal - noisyFunctionReg;
    error2 = computeDerivative(stepsize)*recoveredSignal;
    %Constraints
    endpointConstraint = optimconstr(stepsize - 1);
    for i = 1:stepsize
        if i > 1
            endpointConstraint(i) = recoveredSignal(i) == recoveredSignal(stepsize + i - 1);
        end
    end
    l2problem.Objective = sum(error1.^2) + sum(error2.^2);
    %l2problem.Constraints.cons1 = endpointConstraint;
    l2Solution = solve(l2problem);
end

function l1Solution = SolveL1Problem(noisyFunction, stepsize)
    noisyFunctionReg = zeros(2*stepsize, 1);
    for i = 1:stepsize
        noisyFunctionReg(i) = noisyFunction(i,1);
        noisyFunctionReg(stepsize + i) = noisyFunction(i,2); 
    end
    l1problem = optimproblem('ObjectiveSense', 'min');
    recoveredSignal = optimvar('recoveredSignal', 2*stepsize, 1);
    recoveredSignalDeriv = optimvar('recoveredSignalDeriv', stepsize, 1);
    %Set up terms used in objective
    error1 = recoveredSignal - noisyFunctionReg;
    error2 = computeDerivative(stepsize)*recoveredSignal;
    %Constraints
    endpointConstraint = optimconstr(stepsize - 1);
    for i = 1:stepsize
        if i > 1
            endpointConstraint(i) = recoveredSignal(i) == recoveredSignal(stepsize + i - 1);
        end
    end
    absoluteValueConstraint = optimconstr(2*stepsize);
    for i = 1:2:2*stepsize
        absoluteValueConstraint(i) = -error2((i+1)/2) <= recoveredSignalDeriv((i+1)/2);
        absoluteValueConstraint(i+1) = error2((i+1)/2) <= recoveredSignalDeriv((i+1)/2);
    end
    l1problem.Objective = sum(error1.^2) + sum(recoveredSignalDeriv);
    %l1problem.Constraints.cons1 = endpointConstraint;
    l1problem.Constraints.cons2 = absoluteValueConstraint;
    l1Solution = solve(l1problem);
end

function convertedVector = solutionRegularizer(vector)
    [n, ~] = size(vector);
    convertedVector = zeros(n/2, 2);
    for i = 1:n/2
        convertedVector(i,1) = vector(i);
        convertedVector(i,2) = vector(i + n/2);
    end
end