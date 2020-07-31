%Requires communication and optimization toolbox
opengl software
clear all

%Prompts user to input how many stepsizes should discretize [0,1]
stepsize = input("Enter the amount of steps the function will take: ") + 1;

%Asks user whether they want a random or a sawtooth function
randomOrSawtooth = input("Would you like a sawtooth function (press 1 if yes, otherwise, random function)? ");
if randomOrSawtooth == 1
   
else
signal = buildRandomSignal(stepsize);
end
signalWithError = awgn(signal, 10, 'measured');

%Solve Optimization Problem
problem = optimproblem('ObjectiveSense', 'min');
recoveredSignal = optimvar('recoveredSignal', 1, stepsize);
recoveredSignalDeriv = optimvar('recoveredSignalDeriv', 1, stepsize - 1);
gradient = buildDerivative(stepsize);
error1 = signal - recoveredSignal;
absoluteValueConstraint = optimconstr(2*(stepsize - 1));
recoveredSignalDerivNoAbs = gradient*recoveredSignal';
for i = 1:2:2*(stepsize - 1) - 1
    absoluteValueConstraint(i) = recoveredSignalDerivNoAbs((i+1)/2) <= recoveredSignalDeriv((i+1)/2);
    absoluteValueConstraint(i+1) = -recoveredSignalDerivNoAbs((i+1)/2) <= recoveredSignalDeriv((i+1)/2);
end
problem.Constraints.cons1 = absoluteValueConstraint;
problem.Objective = sum(error1.^2) + sum(recoveredSignalDeriv);
sol = solve(problem);
hold on
plot(0:1/(stepsize-1):1, signal)
plot(0:1/(stepsize-1):1, signalWithError)
plot(0:1/(stepsize-1):1, sol.recoveredSignal)
legend('Original Signal', 'Signal With Error', 'Recovered Signal')
hold off


function signal = buildRandomSignal(stepsize)
    signal = zeros(1, stepsize);
    for i = 1:length(signal)
        signal(i) = rand();
    end
end

function derivativeMatrix = buildDerivative(stepsize)
    derivativeMatrix = zeros(stepsize - 1, stepsize);
    for i = 1:stepsize-1
        derivativeMatrix(i,i) = -1/(stepsize-1);
        derivativeMatrix(i, i+1) = 1/(stepsize - 1);
    end
end