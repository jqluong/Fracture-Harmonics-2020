function plotFaceFunction(V, F, faceFunction)
    %face function is F x 3
    [m,~] = size(faceFunction);
    figure(1)
    hold on
    for i = 1:m
        facePiece = faceFunction(i,:);
        trisurf([1 2 3], V(F(i,:), 1), V(F(i,:), 2), facePiece(:));
    end
    hold off
end