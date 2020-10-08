function Md = face_build_discontinuity_mass(V,F)
    %Want V_d to be (Face1 Vertex1, Face1 Vertex2, Face1 Vertex3, ...)
    Vd = V(F',:);
    %Need to redefine the face matrix.  The first face should be the first
    %three vertices, the second face should be the next three vertices,
    %etc.
    Fd = transpose(reshape(1:3*size(F,1),3,size(F,1)));
    Md = massmatrix(Vd,Fd);
    %There are some Nan's sometimes??
end