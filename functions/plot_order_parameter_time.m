function fg = plot_order_parameter_time(rr)

tt = length(rr);

fg = figure;
plot(rr, '-k', 'LineWidth', 2);
xlim([1 tt]);
xlabel('time');
ylabel('order parameter')
set(gcf,'Position', [3   591   233   151]);
set( gca, 'fontname', 'arial', 'fontsize', 16, 'linewidth', 2 );
box off