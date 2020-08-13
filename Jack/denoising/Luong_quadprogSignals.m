%Requires optimization toolbox
opengl software
clear all

%Prompts user to input how many stepsizes should discretize [0,1]
stepsize = input("Enter the amount of steps the function will take: ");
points = stepsize + 1;
%Only doing square function for now because
signal = 0.5*(buildSquare(stepsize)+1);
%Adds white gaussian noise to the square function
signalWithError = addNoise(signal);
signalWithError = [signalWithError(:,1); signalWithError(:,2)]; 
%Sets up and solves the l^2 optimization problem
%Constraints
alpha = 1;
beta = 1;
gamma = 1;
%Problem Setup
derivativeMatrix = [computeDerivative(stepsize) zeros(stepsize)];
quadraticTerm = beta*(derivativeMatrix' * derivativeMatrix) + alpha*[eye(2*stepsize) zeros(2*stepsize, stepsize); zeros(stepsize,2*stepsize) zeros(stepsize)];
linearTerm = -2*alpha*[signalWithError; zeros(stepsize,1)] + gamma*[zeros(2*stepsize,1); ones(stepsize,1)];
constraintMatrix = [-createEndpointMatrix(stepsize) -eye(stepsize); -createEndpointMatrix(stepsize) eye(stepsize)];
options = optimoptions(@quadprog, 'Algorithm', 'interior-point-convex', 'OptimalityTolerance', 1e-12);
solutionAppend = quadprog(2*quadraticTerm, linearTerm, constraintMatrix, zeros(2*stepsize,1));
%Plotting
solution = [solutionAppend(1:stepsize, 1) solutionAppend(stepsize+1:2*stepsize,1)];
signalWithError = [signalWithError(1:stepsize, 1) signalWithError(stepsize+1:2*stepsize,1)];

hold on
axis tight
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
        p3 = plot([x(i) x(i+1)], [solution(i, 1) solution(i,2)], 'Color', 'g', 'DisplayName', 'Solution', 'LineWidth', 1);
    else
        plot([x(i) x(i+1)], [solution(i, 1) solution(i,2)], 'Color', 'g', 'DisplayName', 'Solution', 'LineWidth', 1)
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