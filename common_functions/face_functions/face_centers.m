function C = face_centers(V, F)
    % this function finds the center coordinates of each face F
    % input is mesh V and F
    % output is #F x dim matrix C
    
    C = zeros(size(F,1), size(V,2));
    for i = 1:length(F)
        v = V(F(i,:),:);
        C(i,:) = mean(v);
    end
end