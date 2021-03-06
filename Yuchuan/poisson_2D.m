%% Solves the 2D Poisson equation: Laplace u = f (over triangular mesh in R^2) subject to Dirichlet boundary condition
%requires gptoolbox
%INPUTS:
    %V = vertices
    %F = faces
    %f = source term of Poisson equation
    %u_0 = boundary condition
    
%OUTPUT:
    %u_int = vector of values of u in the interior of domain

function [u_int] = poisson_2D(V,F,f,u_0)

  % Based on 2.1 of "Algorithms and Interfaces for Real-Time Deformation of 2D and 3D
  % shapes" [Jacobson 2013]

  L = cotmatrix(V,F); %Laplacian operator (sparse)
  M = massmatrix(V,F,'full'); %mass matrix (sparse)
  b = M*f - L*u_0;
  u_int = L\b; 

end
