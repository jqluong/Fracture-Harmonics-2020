function plot_surf(F,V,u)

plot = tsurf(F,[V u],fphong, falpha(1,0));

colormap(cbrewer('RdYlBu',40));
colorbar();

axis equal
grid off;
axis off;
camlight;
add_isolines(plot);

end