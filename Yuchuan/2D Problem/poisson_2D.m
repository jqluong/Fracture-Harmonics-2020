%% Solves the 2D Poisson equation: Laplace u = f (over triangular mesh in R^2) subject to Dirichlet boundary condition
%requires gptoolbox
%INPUTS:
    %V = vertices
    %F = faces
    %f = source term of Poisson equation
    %u_0 = boundary condition
    %bv = boundary vertices
    
%OUTPUT:
    %u_int = vector of values of u in the interior of domain

function [u_sol] = poisson_2D(V,F,f,u_0,bv)

  % Based on 2.1 of "Algorithms and Interfaces for Real-Time Deformation of 2D and 3D
  % shapes" [Jacobson 2013]
  
  
  %if nargin == 5, this is a boundary condition problem

  if nargin == 5
      L = cotmatrix(V,F); %Laplacian operator (sparse)
      M = massmatrix(V,F,'full'); %mass matrix (sparse)

      b_1 = M*f;
      b_1(bv) = [];

      b_2 = (-1)*L*u_0;
      b_2(bv) = [];

      b = b_1 + b_2;


      L(bv,:) = [];
      L(:,bv) = [];

      u_sol = L\b; 

  
  %if nargin == 3, this is a unit norm constraint problem ||u|| = 1    
  elseif nargin == 3
      L = cotmatrix(V,F); %Laplacian operator (sparse)
      M = massmatrix(V,F,'full'); %mass matrix (sparse)
      G  = grad(V,F);
      dblA = doublearea(V,F);
      GMG = G'*repdiag(diag(sparse(dblA)/2),size(V,2))*G;
      %H = inv(diag(sparse(dblA)/2))*GMG;   %how to implement area matrix??
      [eiv,~] = eig(full(GMG));
      u_sol = eiv;
      
  end

end
