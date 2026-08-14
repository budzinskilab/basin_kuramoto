function fg = imagesc_phase(theta, N)

fg = figure;
imagesc(theta');
xlabel('time')
ylabel('nodes');
yticks([1 N]);
colormap('hsv');
clim([-pi pi]);
set(gcf,'Position', [397   565   267   138]);
set( gca, 'fontname', 'arial', 'fontsize', 16, 'linewidth', 2 );