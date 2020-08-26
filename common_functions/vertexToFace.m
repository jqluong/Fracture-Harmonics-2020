function faceFunction = vertexToFace(V, F, vertexFunction)
    faceFunction = zeros(size(F));
    [m,n] = size(F);
    for i = 1:length(vertexFunction)
        for j = 1:m
            for k = 1:n
                if F(j,k) == i
                    faceFunction(j,k) = vertexFunction(i);
                end
            end
        end
    end
end
