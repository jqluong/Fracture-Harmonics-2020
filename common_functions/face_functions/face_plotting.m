function face_plotting(V, F, faceFunction)
    %Face function is 3*|F| x 1
    faceFunction = transpose(reshape(faceFunction, 3, length(faceFunction)/3));
    [m,~] = size(faceFunction);
    figure(1)
    for i = 1:m
        if i == 2
            hold on
        end
        facePiece = faceFunction(i,:);
        trisurf([1 2 3], V(F(i,:), 1), V(F(i,:), 2), facePiece(:));
    end
    colorbar
    hold off
end