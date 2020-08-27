function plotFaceFunction(V, F, faceFunction)
    %face function is F x 3
    [m,~] = size(faceFunction);
    figure(1)
    for i = 1:m
        if i == 1
            hold on
        end
        facePiece = faceFunction(i,:);
        trisurf([1 2 3], V(F(i,:), 1), V(F(i,:), 2), facePiece(:));
    end
    colorbar
    hold off
end