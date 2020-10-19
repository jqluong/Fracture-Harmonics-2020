%Objective: Adds noise to a squarewave on [0,1] with stepsizes determined
%by the user.  Implements l1 and l2 denoising algorithms, and computes the
%gradient of both.  Plots all results.

%Requires communication and optimization toolbox
opengl software
clear all

%Prompts user to input how many stepsizes should discretize [0,1]
stepsize = input("Enter the amount of steps the function will take: ") + 1;

%Asks user whether they want a random or a sawtooth function
randomOrSquare = input("Would you like a square wave (press 1 if yes, otherwise, random function)? ");
if randomOrSquare == 1
    t = 0:1/(stepsize-1):1;
   signal = square(2*pi*t);
else
    signal = buildRandomSignal(stepsize);
end
signalWithError = awgn(signal, 10, 'measured'); %Add noise, 2nd parameter
% is signal to noise ratio

%Set up Optimization Problem
problem1 = optimproblem('ObjectiveSense', 'min');
problem2 = optimproblem('ObjectiveSense', 'min');
recoveredSignal = optimvar('recoveredSignal', 1, stepsize);
recoveredSignalDeriv = optimvar('recoveredSignalDeriv', 1, stepsize - 1);
gradient = buildDerivative(stepsize);
error1 = recoveredSignal - signalWithError;
error2 = gradient*recoveredSignal';
absoluteValueConstraint = optimconstr(2*(stepsize - 1));
recoveredSignalDerivNoAbs = gradient*recoveredSignal';
for i = 1:2:2*(stepsize - 1) - 1
    absoluteValueConstraint(i) = recoveredSignalDerivNoAbs((i+1)/2) <= recoveredSignalDeriv((i+1)/2);
    absoluteValueConstraint(i+1) = -recoveredSignalDerivNoAbs((i+1)/2) <= recoveredSignalDeriv((i+1)/2);
end
problem1.Constraints.cons1 = absoluteValueConstraint;
problem1.Objective = sum(error1.^2) + sum(recoveredSignalDeriv);
problem2.Objective = sum(error1.^2) + sum(error2.^2);

%Solve Optimization Problem
sol1 = solve(problem1);
sol2 = solve(problem2);

%Plot orignal, noisy, and denoised functions
hold on
plot(0:1/(stepsize-1):1, signal)
plot(0:1/(stepsize-1):1, signalWithError)
%plot(0:1/(stepsize-1):1, sol1.recoveredSignal)
%plot(0:1/(stepsize-1):1, sol2.recoveredSignal)
legend('Original Signal', 'Signal With Error', 'Recovered Signal with l^1 term', 'Recovered Signal with l^2 term')
title('Denoising Algorithm with pointwise Functions')
hold off

%Plots gradient of denoised functions
figure
hold on
plot(1:1:stepsize - 1, gradient*sol1.recoveredSignal')
plot(1:1:stepsize - 1, gradient*sol2.recoveredSignal')
legend('Gradient vector of Recovered Signal with l1 term', 'Gradient vector of Recovered Signal with l2 term')
title('Gradient of Denoised Functions')
hold off


%Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Builds a random signal
function signal = buildRandomSignal(stepsize)
    signal = zeros(1, stepsize);
    for i = 1:length(signal)
        signal(i) = rand();
    end
end

%Builds a derivative matrix
function derivativeMatrix = buildDerivative(stepsize)
    derivativeMatrix = zeros(stepsize - 1, stepsize);
    for i = 1:stepsize-1
        derivativeMatrix(i,i) = -1;%*(stepsize-1);
        derivativeMatrix(i, i+1) = 1;%*(stepsize - 1);
    end
end