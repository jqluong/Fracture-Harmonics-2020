function [V,F] = draw_mesh(H,Q,M)

%Input: 
   %H = #H-by-2 list of at least one point in every hole of the domain
        %e.g. if domain is the annulus {z : 1<|z|<2} then H = [0 0] 
   %Q = Quality (as in Triangle documentation)
   %M = MaxArea (as in Triangle documentation)
       
   %Comment (yc): the current code only works for a boundary draw in ONE
   %continuous curve (to be changed)
 
%Output: 
   %F = #F-by-3 faces list
   %V = %V-by-3 vertices list
   
   
%Draw boundary of mesh
P = get_pencil_curve();  %list of boundary points

%Create list of constraint edge indices, E
E1 = (1:1:size(P,1))';
E2 = (2:1:size(P,1)+1)';
E2(size(P,1)) = 1;
E = [E1 E2];

%Triangulate mesh
[V, F] = triangle(P,E,H, 'Quality', Q, 'MaxArea', M, 'Flags','-q1.2-a0.01');



end