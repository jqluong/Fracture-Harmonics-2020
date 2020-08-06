%% problem
freq=1000;
t=0:1/freq:1;

si=zeros(length(t)-1,1);
for i = 1:(length(si)/2)
    si(i) = 1;
end
sr = 0.2*(rand([length(si),1])-0.5);
s=si+sr;
n = length(s);

f = condat_denoise(s, 0.5);

NRMSE = sqrt(mean((si - f).^2))/range(si);
fprintf('normalized root mean square error for l1 signal recovery: %f \n', NRMSE);

hold on
plot(0:1/(n-1):1, s);
plot(0:1/(n-1):1, f);
legend('Noisy Signal', 'Recovered Signal')
hold off

%% 1D Direct TV denoising algorithm from Condat
function x = condat_denoise(y, lambda)
    % https://hal.archives-ouvertes.fr/file/index/docid/675043/filename/condat_killer_tv.pdf
    % params: y = noisy signal, lambda = weighting parameter
    % minimize  (1/2)||x - y||_2^2 + lambda * sum_i |x_{i+1} - x_i|
    % note that sum_i |x_{i+1} - x_i| is l1 norm of 2 point forward finite difference derivative
    n = length(y);
    x = zeros(n,1);
    % part A
    k = 1;
    k0 = 1; % k_0
    km = 1; % k_-
    kp = 1; % k_+
    vmin = y(1) - lambda;
    vmax = y(1) + lambda;
    umin = lambda;
    umax = -lambda;
    
    % part B
    while 1 % elimiate goto in favor of loops
        if k == n
            x(k) = vmin + umin;
            return;
        end
    
        if y(k+1) + umin < vmin - lambda
            x(k0:km) = vmin;
            k = km + 1;
            k0 = k;
            km = k;
            kp = k;
            vmin = y(k);
            vmax = y(k) + 2*lambda;
            umin = lambda;
            umax = -lambda;
        elseif y(k+1) + umax > vmax + lambda
            x(k0:kp) = vmax;
            k = kp + 1;
            k0 = k;
            km = k;
            kp = k;
            vmin = y(k) - 2*lambda;
            vmax = y(k);
            umin = lambda;
            umax = -lambda;
        else
            k = k+1;
            umin = umin + y(k) - vmin;
            umax = umax + y(k) - vmax;
            if umin >= lambda
                vmin = vmin + (umin - lambda)/(k - k0 + 1);
                umin = lambda;
                km = k;
            end
            if umax <= -lambda
                vmax = vmax + (umax + lambda)/(k - k0 + 1);
                umax = -lambda;
                kp = k;
            end
        end
    % part C
        if k >= n
            if umin < 0
                x(k0:km) = vmin;
                k = km + 1;
                k0 = k;
                km = k;
                vmin = y(k);
                umin = lambda;
                umax = y(k) + lambda - vmax;
            elseif umax > 0
                x(k0:km) = vmax;
                k = kp + 1;
                k0 = k;
                kp = k;
                vmax = y(k);
                umin = y(k) - lambda - vmin;
                umax = -lambda;
            else
                x(k0:n) = vmin + umin/(k - k0 + 1);
                return
            end
        end
    end
end