function fg = plot_phase_scatter(N, theta)

theta = angle(exp(1i*theta));

fg = figure;
scatter(1:N, theta, 'k', 'filled');
xlabel('nodes');
ylabel('\Delta \theta');
ylim([-3.5 3.5]);
yticks([-pi, 0, pi]);
yticklabels({'-\pi', '0', '\pi'});
xlim([0 N+1]);
xticks([1 51]);
set(gcf,'Position', [3   591   233   151]);
set( gca, 'fontname', 'arial', 'fontsize', 14, 'linewidth', 2 );
box off