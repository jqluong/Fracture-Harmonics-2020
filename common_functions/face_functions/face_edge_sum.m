function D_sum = face_edge_sum(F)
    %Sums up the discontinuity of a function by summing up the area between
    %the same line segment that share the same face.  Calcultes area by
    %trapezoids.  Used in tandem with discontinuity.m.  To calculate the 
    %discontinuity of a function u defined on the face of a mesh, 
    %call face_edge_sum(F) * discontinuity*(V,F) *u. 
    E = edges(F);
    [m,~] = size(E); %m is number of edges
    D_sum = 1/2*[speye(m,m) speye(m,m)];
end