function dirichletEnergy = computeDirichletEnergy(xaxis, yaxis, f)
    %Computes the Dirichlet Energy of a 2D function f given the xaxis and
    %yaxis the function is defined on. The function f is a matrix based 
    %on xaxis and yaxis, but it gets changed to a column vector later. 
    %Assumes equally spaced partitions.  Assumes you don't have the 
    %triangle mesh yet.
    
    %Area for integration at the end
    cellArea = abs(xaxis(2) - xaxis(1)) * abs(yaxis(2) - yaxis(1));
    
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
    gradient = (G*f(:)) .* (G*f(:));
    dirichletEnergy = (cellArea/2)*sum(gradient, 'all');
end