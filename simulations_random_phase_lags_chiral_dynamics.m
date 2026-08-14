%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                   %%
%% Linear stability and basin size estimatation in Kuramoto networks %%
%% Roberto Budzinski & Kalel Rossi                                   %%
%% August 2026                                                       %%
%%                                                                   %%
%% Simulations networks with random phase-lags and chiral dynamics   %%
%%                                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% setup

clearvars; clc;
addpath( './functions' );

folderPath = 'files_scans/random_phase_lags_data';
if ~isfolder(folderPath)
    mkdir(folderPath);
end

%% parameters

dt = 0.001; T = 25.0; t = 0:dt:T;  %time
N = 51;                            %number of oscillators
epsilon = 1;                       %coupling strength
f_mu = 0;                          %(Hz) natural frequency
omega = ( f_mu*ones(N,1) )*2*pi;   %natural frequency

% adjacency matrix
k = 25;
a = double(ring_graph(N,k));

Ns = 10000;

%% example 1

kk = 10;
M = 40;
[phi, h] = random_circulant_graph_sparse(N, M, kk);
phi = phi.*pi;
for ii = 1:N
    phi(ii,ii) = 0;
end

phi_1 = phi;

%composite matrix
K = (epsilon) .* exp(-1i*phi) .* a;

%eigenspectrum of K
[v,d] = circulant_eigensystem( K ); %analytical eigensystem

ll = repmat(real(diag(d)),3,1);

lambda = nan(N,N);
lambda_max = nan(N,1);

for ii = N+1:2*N
    for jj = 1:N-1
        lambda(ii-N,jj) = (ll(ii-jj) + ll(ii+jj))/2 - ll(ii);
    end
    lambda_max(ii-N) = max(lambda(ii-N,:));
end

lambda_plot = lambda_max;
lambda_plot(floor(N/2):-1:1) = lambda_max(N:-1:floor(N/2)+2);
lambda_plot(floor(N/2)+1) = lambda_max(1);
lambda_plot(floor(N/2)+2:N) = lambda_max(2:floor(N/2)+1);

lambda_plot_1 = lambda_plot;

gamma = diag( real(d) );

for ii = 1:N
    if ( lambda_max(ii) >= 0 )
        gamma(ii) = 0;
    end
end

gamma = power( gamma, N);
gamma = gamma ./ nansum(gamma);

gamma_plot = gamma;
gamma_plot(floor(N/2):-1:1) = gamma(N:-1:floor(N/2)+2);
gamma_plot(floor(N/2)+1) = gamma(1);
gamma_plot(floor(N/2)+2:N) = gamma(2:floor(N/2)+1);
gamma_plot( lambda_plot >= 0 ) = nan;

gamma_plot_1 = gamma_plot;

theta_final  = nan(Ns, N);
rr = nan(Ns, 1);

for ii_s = 1:Ns

    rng(ii_s)
    theta0 = 2*pi*( rand(N,1) - 0.5 );

    %numerical simulation
    [theta_final(ii_s,:) rr(ii_s,:)] = simulate_KM_final_state(a,omega,epsilon,theta0,t,dt,phi);

end

s  = nan(Ns,N);
cc = zeros(N,1);

for jj = 1:N
    for ii = 1:Ns
        s(ii,jj) = order_parameter( theta_final(ii,:) - angle(v(:,jj) )', N );
    end
    cc(jj) = nansum( s(:,jj) > 0.99 );
end

total_q_state_1 = sum(cc);

cc = cc ./ max(cc);
cc = cc ./ nansum(cc);

cc_plot = cc;
cc_plot(floor(N/2):-1:1) = cc(N:-1:floor(N/2)+2);
cc_plot(floor(N/2)+1) = cc(1);
cc_plot(floor(N/2)+2:N) = cc(2:floor(N/2)+1);
cc_plot( lambda_plot >= 0 ) = nan;

cc_plot_1 = cc_plot;

rng(1)
theta0 = 2*pi*( rand(N,1) - 0.5 );

%numerical simulation
[theta_1, ~] = simulate_KM(a,omega,epsilon,theta0,t,dt,phi);

%% example 1

kk = 2135245;
M = 40;
[phi, h] = random_circulant_graph_sparse(N, M, kk);
phi = phi.*pi;
for ii = 1:N
    phi(ii,ii) = 0;
end

phi_2 = phi;

%composite matrix
K = (epsilon) .* exp(-1i*phi) .* a;

%eigenspectrum of K
[v,d] = circulant_eigensystem( K ); %analytical eigensystem

gamma2 = diag(real(d)); %example 1

ll = repmat(real(diag(d)),3,1);
ll_1 = ll;

lambda = nan(N,N);
lambda_max = nan(N,1);

for ii = N+1:2*N
    for jj = 1:N-1
        lambda(ii-N,jj) = (ll(ii-jj) + ll(ii+jj))/2 - ll(ii);
    end
    lambda_max(ii-N) = max(lambda(ii-N,:));
end

lambda_plot = lambda_max;
lambda_plot(floor(N/2):-1:1) = lambda_max(N:-1:floor(N/2)+2);
lambda_plot(floor(N/2)+1) = lambda_max(1);
lambda_plot(floor(N/2)+2:N) = lambda_max(2:floor(N/2)+1);

lambda_plot_2 = lambda_plot;

gamma = diag( real(d) );

for ii = 1:N
    if ( lambda_max(ii) >= 0 )
        gamma(ii) = 0;
    end
end

gamma = power( gamma, N);
gamma = gamma ./ nansum(gamma);

gamma_plot = gamma;
gamma_plot(floor(N/2):-1:1) = gamma(N:-1:floor(N/2)+2);
gamma_plot(floor(N/2)+1) = gamma(1);
gamma_plot(floor(N/2)+2:N) = gamma(2:floor(N/2)+1);
gamma_plot( lambda_plot >= 0  ) = nan;

gamma_plot_2 = gamma_plot;

theta_final  = nan(Ns, N);
rr = nan(Ns, 1);

for ii_s = 1:Ns

    rng(ii_s)
    theta0 = 2*pi*( rand(N,1) - 0.5 );

    %numerical simulation
    [theta_final(ii_s,:) rr(ii_s,:)] = simulate_KM_final_state(a,omega,epsilon,theta0,t,dt,phi);

end

s  = nan(Ns,N);
cc = zeros(N,1);

for jj = 1:N
    for ii = 1:Ns
        s(ii,jj) = order_parameter( theta_final(ii,:) - angle(v(:,jj) )', N );
    end
    cc(jj) = nansum( s(:,jj) > 0.99 );
end

total_q_state_2 = sum(cc);


cc = cc ./ max(cc);
cc = cc ./ nansum(cc);

cc_plot = cc;
cc_plot(floor(N/2):-1:1) = cc(N:-1:floor(N/2)+2);
cc_plot(floor(N/2)+1) = cc(1);
cc_plot(floor(N/2)+2:N) = cc(2:floor(N/2)+1);
cc_plot( lambda_plot >= 0  ) = nan;

cc_plot_2 = cc_plot;

rng(1)
theta0 = 2*pi*( rand(N,1) - 0.5 );

%numerical simulation
[theta_2, ~] = simulate_KM(a,omega,epsilon,theta0,t,dt,phi);

save('files_scans/random_phase_lags_data/figure_3.mat');