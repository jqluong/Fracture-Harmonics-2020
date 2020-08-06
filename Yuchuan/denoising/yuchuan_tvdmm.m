%% Example: TV Denoising using MM
% TV denoising using an algorithm derived using majorization-minimization (MM)
% and solvers for sparse banded systems.
% See accompanying notes.
%
% Ivan Selesnick, selesi@nyu.edu, 2011.
% Revised 2017

%% Start

clear


%% Create data

t = (0:0.01:10)';
clean = square(t); 
noisy = awgn(clean,10,'measured');
                              

figure(1)
clf
subplot(2,1,1)
plot(clean)
title('Clean signal')

subplot(2,1,2)
plot(noisy)
title('Noisy signal')




%% TV Denoising
% Run TV denoising algorithm (MM algorithm)

lam = 1.5;                         % lam: regularization parameter
Nit = 10;                          % Nit: number of iterations
[x, cost] = tvd_mm(noisy, lam, Nit);   % Run MM TV denoising algorithm



figure(1)
clf
subplot(3,1,1)
plot(clean)
title('Clean signal')

subplot(3,1,2)
plot(noisy)
title('Noisy signal')

subplot(3,1,3)
plot(x)
title('TV denoising')

