%% Create mesh
[x, y] = meshgrid(1:10, 1:10);
points = [x(:), y(:)];
V = [points zeros(size(points,1),1)];
F = delaunay(points);

%% Eigen-problem
M_inv = inv(face_area_matrix(V,F));
GMG = face_GMG(V,F);
H = full((1/2) * M_inv * GMG);
[eiv,eit] = eig(H);


%% Plotting

face_plotting(V,F,eiv(:,164)); %plot 1,2,164 and 325






