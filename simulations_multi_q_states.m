%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                   %%
%% Linear stability and basin size estimatation in Kuramoto networks %%
%% Roberto Budzinski & Kalel Rossi                                   %%
%% August 2026                                                       %%
%%                                                                   %%
%% Simulations of Kuramoto networks for multi-q-states               %%
%%                                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% setup

clearvars; clc;
addpath( './functions' );

folderPath = 'files_scans/multi_q_states_data';
if ~isfolder(folderPath)
    mkdir(folderPath);
end

%% parameters

dt = 0.001; T = 50.0; t = 0:dt:T;   %time
N = 51;                             %number of oscillators
epsilon = 10;                       %coupling strength
f_mu = 0;                           %(Hz) natural frequency
omega = ( f_mu*ones(N,1) )*2*pi;    %natural frequency
Ttrans = 10000;                     %transient discarded before order-parameter stats

%% example 1: k = 4

ii_k  = 4;

% adjacency matrix
k = ii_k;
a = double(ring_graph(N,k));

%phase-lag
phi = heterogeneous_phase_lag(a,pi,k);

%composite matrix
K = (epsilon) .* exp(-1i*phi) .* a;

%eigenspectrum of K
[v,d] = circulant_eigensystem( K );

%numerical simulations

Ns = 1000; %number of realizations

theta_final = nan(Ns, N);
rr = nan(Ns, length(t));

for ii_s = 1:Ns

    %random initial conditions
    rng(ii_s)
    theta0 = 2*pi*( rand(N,1) - 0.5 );

    [theta_aux, rr(ii_s,:)] = simulate_KM(a, omega, epsilon, theta0, t, dt, phi);
    theta_final(ii_s,:) = theta_aux(end,:);
    clear theta_aux;

end

%order-parameter summary stats over the post-transient window (avoids saving the full rr trace)
order_mean = nan(Ns,1);
order_std = nan(Ns,1);

for ii = 1:Ns
    order_mean(ii) = nanmean(rr(ii,Ttrans:end));
    order_std(ii) = nanstd(rr(ii,Ttrans:end));
end

clear rr;

str = sprintf('files_scans/multi_q_states_data/N_%d_k_%d.mat', N, k);
save(str);

%% example 2: k = 6

ii_k  = 6;

% adjacency matrix
k = ii_k;
a = double(ring_graph(N,k));

%phase-lag
phi = heterogeneous_phase_lag(a,pi,k);

%composite matrix
K = (epsilon) .* exp(-1i*phi) .* a;

%eigenspectrum of K
[v,d] = circulant_eigensystem( K );

%numerical simulations

Ns = 1000; %number of realizations

theta_final = nan(Ns, N);
rr = nan(Ns, length(t));

for ii_s = 1:Ns

    %random initial conditions
    rng(ii_s)
    theta0 = 2*pi*( rand(N,1) - 0.5 );

    [theta_aux, rr(ii_s,:)] = simulate_KM(a, omega, epsilon, theta0, t, dt, phi);
    theta_final(ii_s,:) = theta_aux(end,:);
    clear theta_aux;

end

%order-parameter summary stats over the post-transient window (avoids saving the full rr trace)
order_mean = nan(Ns,1);
order_std = nan(Ns,1);

for ii = 1:Ns
    order_mean(ii) = nanmean(rr(ii,Ttrans:end));
    order_std(ii) = nanstd(rr(ii,Ttrans:end));
end

clear rr;

str = sprintf('files_scans/multi_q_states_data/N_%d_k_%d.mat', N, k);
save(str);
