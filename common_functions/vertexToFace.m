function faceFunction = vertexToFace(F, vertexFunction)
    [m,n] = size(F);
    faceFunction = zeros(m*n,1);
    for i = 1:length(vertexFunction)
        for j = 1:m
            for k = 1:n
                if F(j,k) == i
                    faceFunction(3*(j-1) + k) = vertexFunction(i);
                end
            end
        end
    end
end
