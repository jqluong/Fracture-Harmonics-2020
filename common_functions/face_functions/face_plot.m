function face_plot(V, F, u)
    % plot face valued function
    %Face function is 3*|F| x 1
    
    u = transpose(reshape(u, 3, length(u)/3));
    [m,~] = size(u);
    figure(1)
    for i = 1:m
        if i == 2
            hold on
        end
        facePiece = u(i,:);
        trisurf([1 2 3], V(F(i,:), 1), V(F(i,:), 2), facePiece(:));
    end
    colormap(cbrewer('RdYlBu',40));
    colorbar
    hold off
end