%Objective: Generates a square signal, adds noise to it, and tries to
%recover the original signal via the minimization of ||f-g||^2 + ||Gf||^2 +
%||Df||_1.  Treats functions as collection of line segments.  G is the
%gradient matrix and D is the discontinuity matrix.
%Functions are defined as a 2n x 1 vector where for each interval i out of n, 
%the startpoint is given at location i and the endpoint is given at location n
%+ i.  They are also defined as n x 2 vector where the first column
%contains startpoints and the second column contains endpoints.

%Requires optimization toolbox
opengl software
clear all

%Prompts user to input how many segmentsshould discretize [0,1]
segmentsnumber = input("Enter the amount of steps the function will take: ");
points = segmentsnumber + 1;
%Generates square function
signal = 0.5*(buildSquare(segmentsnumber)+1);
%Adds white gaussian noise to the square function
signalWithError = addNoise(signal);
signalWithError = [signalWithError(:,1); signalWithError(:,2)]; 
%Sets up and solves the optimization problem via quadprog.
%Constraints
alpha = 1;
beta = 1;
gamma = 1;
%Problem Setup
derivativeMatrix = [computeDerivative(segmentsnumber) sparse(segmentsnumber, segmentsnumber)];
quadraticTerm = beta*(derivativeMatrix' * derivativeMatrix) + alpha*[speye(2*segmentsnumber,2*segmentsnumber) sparse(2*segmentsnumber, segmentsnumber); sparse(segmentsnumber,2*segmentsnumber) sparse(segmentsnumber,segmentsnumber)];
linearTerm = -2*alpha*[signalWithError; zeros(segmentsnumber,1)] + gamma*[zeros(2*segmentsnumber,1); ones(segmentsnumber,1)];
constraintMatrix = [-createEndpointMatrix(segmentsnumber) -speye(segmentsnumber); createEndpointMatrix(segmentsnumber) -speye(segmentsnumber)];
options = optimoptions(@quadprog, 'Algorithm', 'interior-point-convex', 'OptimalityTolerance', 1e-12);
solutionAppend = quadprog(2*quadraticTerm, linearTerm, constraintMatrix, zeros(2*segmentsnumber,1));
%Plotting
solution = [solutionAppend(1:segmentsnumber, 1) solutionAppend(segmentsnumber+1:2*segmentsnumber,1)];
signalWithError = [signalWithError(1:segmentsnumber, 1) signalWithError(segmentsnumber+1:2*segmentsnumber,1)];

hold on
axis tight
x = 0:1/segmentsnumber:1;
for i = 1:segmentsnumber
    if i == 1
        p1 = plot([x(i) x(i+1)], [signal(i, 1) signal(i,2)], 'Color', 'b', 'DisplayName', 'Original Signal', 'LineWidth', 1);
    else
        plot([x(i) x(i+1)], [signal(i, 1) signal(i,2)], 'Color', 'b', 'DisplayName', 'Original Signal', 'LineWidth', 1)
    end
end
%for i = 1:segmentsnumber
%    if i == 1
%        p2 = plot([x(i) x(i+1)], [signalWithError(i, 1) signalWithError(i,2)], 'Color', 'r', 'DisplayName', 'Noisy Signal', 'LineWidth', 1);
%    else
%        plot([x(i) x(i+1)], [signalWithError(i, 1) signalWithError(i,2)], 'Color', 'r', 'DisplayName', 'Noisy Signal', 'LineWidth', 1)
%    end
%end
for i = 1:segmentsnumber
    if i == 1
        p3 = plot([x(i) x(i+1)], [solution(i, 1) solution(i,2)], 'Color', 'g', 'DisplayName', 'Solution', 'LineWidth', 1);
    else
        plot([x(i) x(i+1)], [solution(i, 1) solution(i,2)], 'Color', 'g', 'DisplayName', 'Solution', 'LineWidth', 1)
    end
end
%for i = 1:segmentsnumber
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

%Builds an n x 2n derivative matrix for [startpoint; endpoint] functions
%using sparse matrices
function gradientMatrix = computeDerivative(segmentsnumber)
    gradientMatrix = [-speye(segmentsnumber) speye(segmentsnumber)];
end

%Builds a n x 2n matrix to evaluate endpoint differences using sparpse
%matrices
function endpointMatrix = createEndpointMatrix(segmentsnumber)
    
    endpointMatrix = [sparse(1, segmentsnumber*2); sparse(segmentsnumber-1,1) -speye(segmentsnumber-1) speye(segmentsnumber-1) sparse(segmentsnumber-1,1)];
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