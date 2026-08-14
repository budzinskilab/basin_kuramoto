%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                   %%
%% Linear stability and basin size estimatation in Kuramoto networks %%
%% Roberto Budzinski & Kalel Rossi                                   %%
%% August 2026                                                       %%
%%                                                                   %%
%% Figure 1 — stability analysis                                     %%
%%                                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% setup

clearvars; clc;
addpath( './functions' );

%% parameters

dt = 0.001; T = 2.0; t = 0:dt:T;    %time
N = 51;                             %number of oscillators
epsilon = 1;                        %coupling strength
f_mu = 0;                           %(Hz) natural frequency
omega = ( f_mu*ones( N,1 ) )*2*pi;  %natural frequency


%% example stability with and without phase-lags

%initial condition for simulation
rng(0);
theta0 = 2*pi*( rand( N,1 ) - 0.5 );

% adjacency matrix
k = 20;
A = double( ring_graph( N, k ) );

% no phase-lag
phi = zeros(N);

lambda_no_pl = linear_stability( A, epsilon, phi ); %stability

theta_no_pl = simulate_KM( A, omega, epsilon, theta0, t, dt, phi ); %simulation

% phase-lags
phi = heterogeneous_phase_lag( A, pi, k );

lambda_pl = linear_stability( A, epsilon, phi ); %stability

theta_pl = simulate_KM( A, omega, epsilon, theta0, t, dt, phi ); %simulation

%plots

fg1 = figure; hold on
stem( lambda_no_pl, 'filled');
stem( lambda_pl, 'filled');
ylabel('\lambda');
xlabel('q-state');
xticks([6 16 26 36 46])
xticklabels({'-20', '-10', '0', '10', '20'});
lb = legend({'no phase-lag', 'phase-lag'});
lb.EdgeColor = 'white';
lb.Location = 'southeast';
lb.FontSize = 10;
set(gcf,'Position', [ 273   760   391   213]);
set( gca, 'fontname', 'arial', 'fontsize', 13, 'linewidth', 2 );

fg2 = figure;
imagesc(t, 1:N, theta_no_pl');
xlabel('time')
ylabel('nodes')
yticks([1 51]);
xticks([0 1]);
xlim([0 1.500])
colormap('hsv');
clim([-pi pi]);
set(gcf,'Position', [129   565   267   138]);
set( gca, 'fontname', 'arial', 'fontsize', 12, 'linewidth', 2 );

fg3 = figure;
imagesc(t, 1:N, theta_pl');
xlabel('time')
ylabel('nodes');
yticks([1 51]);
xticks([0 1]);
xlim([0 1.500])
colormap('hsv');
clim([-pi pi]);
set(gcf,'Position', [397   565   267   138]);
set( gca, 'fontname', 'arial', 'fontsize', 12, 'linewidth', 2 );

%% stability phase-lags as a function of k

lambda = nan(floor(N/2), N);

for ii_k = 1:floor(N/2)

    k = ii_k;

    % adjacency matrix
    A = double( ring_graph( N, k ) );

    %phase-lags
    phi = heterogeneous_phase_lag( A, pi, k );

    lambda(ii_k,:) = linear_stability( A, epsilon, phi );

end

%plot

fg4 = figure;
imagesc(lambda'<0);
ylim([1 26.5]);
ylabel('|q|');
xlabel('k');
yticks([1 6 11 16 21 26]);
yticklabels({'25', '20', '15', '10', '5', '0'});
xticks([1 5 10 15 20 25]);
colormap(flipud(bone));
clim([0 1.25]);
set(gcf,'Position', [854   436   640   355]);
grid on
set( gca, 'fontname', 'arial', 'fontsize', 15, 'linewidth', 2 );