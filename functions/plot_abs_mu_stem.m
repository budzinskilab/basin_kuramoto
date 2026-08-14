function fg = plot_abs_mu_stem(N, mu)

mu_plot = mu;
mu_plot(floor(N/2):-1:1) = mu(N:-1:floor(N/2)+2);
mu_plot(floor(N/2)+1) = mu(1);
mu_plot(floor(N/2)+2:N) = mu(2:floor(N/2)+1);

fg = figure;
stem(abs(mu_plot), 'LineWidth', 1, 'markersize', 3, 'Color', [0 0 0],  'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
ylabel('|\mu|');
xlabel('q-states');
xlim([0 N+1]);
xticks([6 16 26 36 46])
xticklabels({'-20', '-10', '0', '10', '20'});
set(gcf,'Position', [237   591   233   151]);
set( gca, 'fontname', 'arial', 'fontsize', 14, 'linewidth', 2 );
box off

