%% Computes a diagonal matrix with edge lengths on the diagonal

%INPUTS
    %V = vertices
    %F = faces
%OUTPUTS
    %EL = |E|-by-|E| matrix
    
function EL = face_edge_lengths_matrix(V,F)
    
    E = edges(F);
    d = zeros(size(E,1),1);
    
    for i = 1:size(E,1)
        v1 = V(E(i,1),:);
        v2 = V(E(i,2),:);
        d(i) = sqrt((v1(1)-v2(1))^2 + (v1(2)-v2(2))^2 + (v1(3)-v2(3))^2);
    end
    
    %EL = diag(d);
    EL = sparse(1:length(d), 1:length(d), d);
end