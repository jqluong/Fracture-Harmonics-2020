function dirichletEnergy = computeDirichletEnergy(xaxis, yaxis, f)
    %Computes the Dirichlet Energy of a 2D function f given the xaxis and
    %yaxis the function is defined on. The function f is a matrix based 
    %on xaxis and yaxis, but it gets changed to a column vector later. 
    %Assumes you don't have the triangle mesh yet.
    
    %Need to do this step and idk why
    [xaxis, yaxis] = meshgrid(xaxis, yaxis);
    %gptoolbox grad function requires vertex matrix to be n x 2
    xaxis = xaxis(:);
    yaxis = yaxis(:);
    %Creating the face matrix
    T = delaunay(xaxis,yaxis);
    %Constructs gradient matrix
    G = grad([xaxis yaxis], T); 
    %Compute dirichlet Energy
    gradient = (G(1:length(G)/2, :)*f(:)).^2 + (G(length(G)/2+1:length(G), :)*f(:)).^2;
    cellArea = doublearea([xaxis yaxis], T)/2;
    %Weigh each gradient term with the cellarea
    dirichletEnergy = sum(cellArea'*gradient, 'all');
end