%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                   %%
%% Linear stability and basin size estimatation in Kuramoto networks %%
%% Roberto Budzinski & Kalel Rossi                                   %%
%% August 2026                                                       %%
%%                                                                   %%
%% Figure 3 — networks with random phase-lags and chiral dynamics    %%
%%                                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('files_scans/random_phase_lags_data/figure_3.mat');

fg1 = figure; hold on
stem( lambda_plot_1, 'filled');
stem( lambda_plot_2, 'filled');
ylabel('\lambda');
xlabel('q-state');
xticks([6 16 26 36 46])
xticklabels({'-20', '-10', '0', '10', '20'});
ylim([-8 21]);
lb = legend({'matrix 1', 'matrix 2'});
lb.EdgeColor = 'white';
lb.Location = 'southeast';
lb.FontSize = 10;
set(gcf,'Position', [ 471   579   643   213]);
set( gca, 'fontname', 'arial', 'fontsize', 13, 'linewidth', 2 );

fg2 = figure; hold on
stem(cc_plot_1, 'filled', 'MarkerSize', 10, 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5]);
stem(gamma_plot_1, 'filled', 'Color', [0 0 0],  'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
ylabel('relative basin size')
xlabel('q-state');
lb = legend('numerical', 'prediction');
lb.EdgeColor = 'white';
lb.Location = 'north';
ylim([-0.05 1.05]);
%yticks([0 0.25 0.5]);
xlim([0 52])
xticks([6 16 26 36 46])
xticklabels({'-20', '-10', '0', '10', '20'});
set(gcf,'Position', [484   333   293   189]);
set( gca, 'fontname', 'arial', 'fontsize', 13, 'linewidth', 2 );
box off

fg3 = figure; hold on
stem(cc_plot_2, 'filled', 'MarkerSize', 10, 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5]);
stem(gamma_plot_2, 'filled', 'Color', [0 0 0],  'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
ylabel('relative basin size')
xlabel('q-state');
lb = legend('numerical', 'prediction');
lb.EdgeColor = 'white';
lb.Location = 'north';
ylim([-0.05 1.05]);
%yticks([0 0.25 0.5]);
xlim([0 52])
xticks([6 16 26 36 46])
xticklabels({'-20', '-10', '0', '10', '20'});
set(gcf,'Position', [ 796   333   293   189]);
set( gca, 'fontname', 'arial', 'fontsize', 13, 'linewidth', 2 );
box off

fg4 = figure;
imagesc(phi_1)
xlabel('node i');
xticks([1 25 51]);
yticks([1 25 51]);
ylabel('node j');
colormap('parula')
clim([0 pi]);
set(gcf,'Position', [516   851   225   178]);
set( gca, 'fontname', 'arial', 'fontsize', 13, 'linewidth', 2 );

fg5 = figure;
imagesc(phi_2)
xlabel('node i');
ylabel('node j');
xticks([1 25 51]);
yticks([1 25 51]);
colormap('parula')
clim([0 pi]);
set(gcf,'Position', [808   851   225   178]);
set( gca, 'fontname', 'arial', 'fontsize', 13, 'linewidth', 2 );

fg6 = figure;
imagesc(t, 1:N,theta_1')
xlabel('time')
ylabel('nodes')
yticks([1 51]);
xlim([0 5])
colormap('hsv');
clim([-pi pi]);
set(gcf,'Position', [494   108   293   163]);
set( gca, 'fontname', 'arial', 'fontsize', 12, 'linewidth', 2 );

fg7 = figure;
imagesc(t, 1:N,theta_2')
xlabel('time')
ylabel('nodes')
yticks([1 51]);
xlim([1 5])
colormap('hsv');
clim([-pi pi]);
set(gcf,'Position', [799   113   293   163]);
set( gca, 'fontname', 'arial', 'fontsize', 12, 'linewidth', 2 );