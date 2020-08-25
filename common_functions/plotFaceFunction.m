function plotFaceFunction(V, F, faceFunction)
    [m,~] = size(faceFunction);
    for i = 1:m
        facePiece = faceFunction(i,:);
        trisurf(F(i,:), V(F(i,:), 1), V(F(i,:), 2), facePiece(:));
    end
end