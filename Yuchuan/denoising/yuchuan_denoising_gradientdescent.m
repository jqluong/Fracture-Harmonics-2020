clear

t = (0:0.01:10)';

c = 10*square(t); 
clean = [c;c];
noisy = awgn(c,10,'measured');
g = [noisy;noisy];

f_old = ones(length(g),1); %initial guess

N = 10000; %# of iterations
k = 0.2;


for i=1:N
    d2f = [0;diff(f_old,2);0];
    grad_f = k*(f_old - g - d2f);
    f_new = f_old - grad_f;
    f_old = f_new;
    

end

figure(1)
clf
subplot(3,1,1)
plot(clean)
title('Clean signal')
%axis([0 10 -20 20])

subplot(3,1,2)
plot(g)
title('Noisy signal')
%axis([0 10 -20 20])

subplot(3,1,3)
plot(f_old)
title('Denoised')
%axis([0 10 -20 20])






