function x = tvdmm_segments(y, lam, Nit, h)
% [x, cost] = tvd_mm(y, lam, Nit)
% Total variation denoising using majorization-minimization
% and banded linear systems.
%
% INPUT
%   y - noisy signal
%   lam - regularization parameter
%   Nit - number of iterations
%
% OUTPUT
%   x - denoised signal
%   cost - cost function history
%
% Reference
% 'On total-variation denoising: A new majorization-minimization
% algorithm and an experimental comparison with wavalet denoising.'
% M. Figueiredo, J. Bioucas-Dias, J. P. Oliveira, and R. D. Nowak.
% Proc. IEEE Int. Conf. Image Processing, 2006.

% Ivan Selesnick, selesi@nyu.edu, 2011
% Revised 2017

y = y(:);                                              % Make column vector                                % Cost function history
N = length(y);

P = (1/h)*[-eye((N-2)/2) zeros((N-2)/2,2) eye((N-2)/2)];  %N is always even %P is (N-2)/2 by N
PPT = P * P';

x = y;                                                 % Initialization
Px = P*x;
Py = P*y;

for k = 1:Nit
    F = sparse(1:(N-2)/2, 1:(N-2)/2, abs(Px)/lam) + PPT;       % F : Sparse banded matrix
    x = y - P'*(F\Py);                                 % Solve banded linear system
    Px = P*x;
end
