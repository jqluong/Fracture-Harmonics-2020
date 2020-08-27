function G = face_grad(V, F)
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
    % 
    % to get output ordered by faces rather than by vertices, we construct
    % a permutation matrix and multiply G by it.
    % 
    % so the output matrix G, operates on a function u to give:
    % G*u = [du_x(F1) du_y(F1) ... ]'
    
    nF = length(F);
    d = size(V, 2);
    
    V_n = zeros(3*nF, d);
    F_n = reshape(1:3*nF, 3, [])';

    for i = 1:size(F,1)
        j = 3*i - 2;
        V_n(j: j+2,:) = V(F(i,:),:);
    end
    
    G = grad(V_n, F_n);
    
    P = zeros(nF*d, 2);
    p = [1:d:nF*d; 1:nF]';
    P(1:nF,:) = p;
    for i = 1:d-1
        P(i*nF+1:(i+1)*nF,:) =  p + i*[1, nF];
    end
    
    P = sparse(P(:,1), P(:,2), ones(length(P),1));
    G = P*G;
end