function face_plot(V,F,u)
    % plot face valued function
    % this is slow becuase it has to duplicate elments of V as like the
    % gradient
    
    F_n = reshape(1:numel(F),3, [])';

    V_n = zeros(3*size(F,1), size(V,2));
    for i = 1:length(F)
        j = 3*i - 2;
        V_n(j: j+2,:) = V(F(i,:)',:);
    end
    
    tsurf(F_n, [V_n u], fphong, falpha(1,0));
    
    colormap(cbrewer('RdYlBu',40));
    colorbar();
    
    axis equal
    grid off;
    axis off;
end