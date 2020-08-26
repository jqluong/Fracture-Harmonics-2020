function G = grad_face(V, F)
    % find gradient matrix for face defined functions.
    % input: 
    %        vertex matrix V, face matrix F
    % output:
    %        G = #faces*dim by #faces*3 Gradient operator
    %
    % for representations of functions as:
    % u = [ u(F1_v1) u(F1_v2) u(F1_v3) u(F2_v1) u(F2_v2) u(F2_v3) ... ]'
    %
    % first, create new V (V_n) that has duplicate coordinates to match
    % representation of u. Then, have F_n that is face matrix indexing this
    % new V_n. And then call grad(V_n, F_n)
    
    V_n = zeros(3*size(F,1), size(V,2));
    F_n = reshape(1:numel(F),3, [])';

    for i = 1:size(F,1)
        j = 3*i - 2;
        V_n(j: j+2,:) = V(F(i,:)',:);
    end
    G = grad(V_n, F_n);
end