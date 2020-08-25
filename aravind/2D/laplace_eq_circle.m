th = [0:2*pi()/1024:2*pi()]';
x = cos(th);
y = sin(th);

[V, F] = triangle([x,y], 'Quality', 30, 'MaxArea', 0.001);

% Boundary conditions: on circle. if y < 0, u = x^2, if y > 0, u = y^2
BP1 = V(V(:,1).^2 + V(:,2).^2 == 1 & V(:,2)< 0, :);
BP2 = V(V(:,1).^2 + V(:,2).^2 == 1 & V(:,2)> 0, :);
BP = [BP1; BP2];

g = [BP1(:,1).^2; BP2(:,2).^2];

u = laplace_eq_2D(V, F, [BP, g]);

trisurf(F, V(:,1), V(:,2), zeros(size(V(:,1))), u, 'LineStyle','none');
colorbar();
title('\nabla^2u = 0, with u = x^2 at y < 0 and u = y^2 at y > 0 on unit circle ');
xlabel('x \in [-1,1]');
ylabel('y \in [-1,1]');