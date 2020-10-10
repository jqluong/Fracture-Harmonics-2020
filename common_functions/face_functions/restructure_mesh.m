function [V_n, F_n] = restructure_mesh(V, F)
    % restructure the mesh to have #V_n = 3*#F, so that you can easily use
    % the new mesh in prexisting functions
    nF = length(F);
    d = size(V, 2);
    
    V_n = zeros(3*nF, d);
    F_n = reshape(1:3*nF, 3, [])';

    for i = 1:size(F,1)
        j = 3*i - 2;
        V_n(j: j+2,:) = V(F(i,:),:);
    end
end