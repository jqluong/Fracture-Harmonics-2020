[V, F] = mesh_square(10);
B = unique(reshape(outline(F),[],1));
B1 = B(V(B,1) < 0.1);
B2 = B(V(B,1) > 0.95);
B = [B1; B2];
g1 = zeros(length(B1),1);
g2 = ones(length(B2),1);
g = [g1;g2];