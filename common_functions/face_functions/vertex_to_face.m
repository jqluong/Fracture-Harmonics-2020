function f_face = vertex_to_face(F, f_vertex)
    [m,n] = size(F);
    f_face = zeros(m*n,1);
    for i = 1:length(f_vertex)
        for j = 1:m
            for k = 1:n
                if F(j,k) == i
                    f_face(3*(j-1) + k) = f_vertex(i);
                end
            end
        end
    end
end
