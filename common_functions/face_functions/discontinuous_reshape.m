function [Vd,Fd] = discontinuous_reshape(V,F)
% transforms mesh matrices from continuous mesh space to discontinuous mesh
% space
% V: continuous mesh vertices matrix
% F: continuous mesh faces matrix

% transform vertex matrix into discontinuous mesh space
Vd = V(transpose(F),:);

[m,~] = size(F);

% transform faces matrix into discontinuous mesh space
Fd = transpose(reshape(1:3*m,3,m));

end