%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                   %%
%% Linear stability and basin size estimatation in Kuramoto networks %%
%% Roberto Budzinski, Kalel Rossi                                    %%
%% August 2026                                                       %%
%%                                                                   %%
%% Figure 2 — analysis of multi-q-states                             %%
%%                                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% setup

clearvars; clc;
addpath( './functions' );

%% example 1: k = 4

ii_jj = 4;

str = sprintf('files_scans/multi_q_states_data/N_51_k_%d.mat', ii_jj);
load(str);

fg = figure;
h = histogram(order_mean, 100);
xx = h.BinEdges;
yy = h.Values;
close(fg);

fg = figure;
bar(xx(1:end-1), yy, 5);
set(gcf,'Position', [761   503   520   200]);
set( gca, 'fontname', 'arial', 'fontsize', 15, 'linewidth', 2 );
xlabel('mean order parameter')
xlim([-0.0025 0.031])
ylabel('count');
box off

ii = 4;

theta_diff = nan(N,1);

for jj = 2:N
    theta_diff(jj) = theta_final(ii,jj) - theta_final(ii,jj-1);
end
theta_diff(1) = theta_final(ii,1) - theta_final(ii,N);

plot_phase_scatter(N, theta_diff);
set(gcf,'Position', [3   320   253   140]);

mu = nan(1,N);
for jj = 1:N
    mu(jj) = exp(1i*theta_final(ii,:))*v(:,jj);
end

mu = power(abs(mu),2);
mu = mu./N;

plot_abs_mu_stem(N, mu);
set(gcf,'Position', [ 5   124   253   140]);

ii = 32;

theta_diff = nan(N,1);

for jj = 2:N
    theta_diff(jj) = theta_final(ii,jj) - theta_final(ii,jj-1);
end
theta_diff(1) = theta_final(ii,1) - theta_final(ii,N);

plot_phase_scatter(N, theta_diff);
set(gcf,'Position', [ 257   320   253   140]);

mu = nan(1,N);
for jj = 1:N
    mu(jj) = exp(1i*theta_final(ii,:))*v(:,jj);
end

mu = power(abs(mu),2);
mu = mu./N;

plot_abs_mu_stem(N, mu)
set(gcf,'Position', [257   123   253   140]);

ii = 75;

theta_diff = nan(N,1);

for jj = 2:N
    theta_diff(jj) = theta_final(ii,jj) - theta_final(ii,jj-1);
end
theta_diff(1) = theta_final(ii,1) - theta_final(ii,N);

plot_phase_scatter(N, theta_diff);
set(gcf,'Position', [507   317   253   140]);

mu = nan(1,N);
for jj = 1:N
    mu(jj) = exp(1i*theta_final(ii,:))*v(:,jj);
end

mu = power(abs(mu),2);
mu = mu./N;

plot_abs_mu_stem(N, mu);
set(gcf,'Position', [507   120   253   140]);

%% example 2: k = 6

ii_jj = 6;

str = sprintf('files_scans/multi_q_states_data/N_51_k_%d.mat', ii_jj);
load(str);

fg = figure;
h = histogram(order_mean, 98);
xx = h.BinEdges;
yy = h.Values;
close(fg);

fg = figure;
bar(xx(1:end-1), yy, 5);
set(gcf,'Position', [761   503   520   200]);
set( gca, 'fontname', 'arial', 'fontsize', 15, 'linewidth', 2 );
xlabel('mean order parameter')
xlim([-0.0025 0.04])
ylim([-10 700])
ylabel('count');
box off

a = double(ring_graph(N, ii_jj));
phi = heterogeneous_phase_lag(a,pi,ii_jj);

t = 0:dt:T;

ii = 667;

rng(ii)
theta0 = 2*pi*( rand(N,1) - 0.5 );

[theta, rr]= simulate_KM(a,omega,epsilon,theta0,t,dt,phi);

plot_order_parameter_time(rr);
ylim([-0.0025 0.1])
set(gcf,'Position', [806   534   454   221]);

plot_order_parameter_time(rr(end-5000:end));
set(gcf,'Position', [806   534   454   221]);
xlabel(''); ylabel(''); xticks([]); yticks([])

imagesc_phase(theta, N);
set(gcf,'Position', [806   251   647   226]);

theta_diff = nan(length(t), N);

for tt = 1:length(t)
    for jj = 2:N
        theta_diff(tt,jj) = theta(tt,jj) - theta(tt,jj-1);
    end
    theta_diff(tt,1) = theta(tt,1) - theta(tt,N);
end

%evaluating similarity S
ss = nan(length(t), N);

for ii = 1:N
    for tt = 1:length(t)
        ss(tt,ii) = order_parameter((theta(tt,:) - angle(v(:,ii))'),N);
    end
end

mu = nan(length(t), N);
for tt = 1:length(t)
    for jj = 1:N
        mu(tt,jj) = exp(1i*theta(tt,:))*v(:,jj);
    end
end

mu = power(abs(mu),2);
mu = mu./N;

fg = figure; hold on
plot(mu(:,4), 'LineWidth', 2);
plot(mu(:,5), 'LineWidth', 2);
plot(mu(:,6), 'LineWidth', 2);
plot(mu(:,7), 'LineWidth', 2);
plot(mu(:,46), 'LineWidth', 2);
plot(mu(:,47), 'LineWidth', 2);
plot(mu(:,48), 'LineWidth', 2);
plot(mu(:,49), 'LineWidth', 2);
lb = legend('q = 3', 'q = 4', 'q = 5', 'q = 6', 'q = -3', 'q = -4', 'q = -5', 'q = -6');
lb.EdgeColor = 'white';
lb.Location ='northeast';
lb.FontSize = 10;
lb.NumColumns = 3;
xlabel('time');
ylabel('\mu');
ylim([-0.05 0.8]);
xlim([0 length(t)]);
set(gcf,'Position', [806   251   647   226]);
set( gca, 'fontname', 'arial', 'fontsize', 15, 'linewidth', 2 );
