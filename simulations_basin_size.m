%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                   %%
%% Linear stability and basin size estimatation in Kuramoto networks %%
%% Roberto Budzinski & Kalel Rossi                                   %%
%% August 2026                                                       %%
%%                                                                   %%
%% Calculation eigenvalues of K and simulations of Kuramoto networks %%
%%                                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% setup

clearvars -except ii_k; clc;
addpath( './functions' );

folderPath = 'files_scans/basin_size_data';
if ~isfolder(folderPath)
    mkdir(folderPath);
end

%% parameters

dt = 0.001; T = 20.0; t = 0:dt:T;   %time
N = 51;                             %number of oscillators
epsilon = 10;                       %coupling strength
f_mu = 0;                           %(Hz) natural frequency
omega = ( f_mu*ones(N,1) )*2*pi;    %natural frequency

vector_k = linspace(1,floor(N/2), floor(N/2));

% adjacency matrix
k = ii_k;
a = double(ring_graph(N,k));

%phase-lag
phi = heterogeneous_phase_lag(a,pi,k);

%composite matrix
K = (epsilon) .* exp(-1i*phi) .* a;

%eigenspectrum of K
[v,d] = circulant_eigensystem( K );

gamma = diag(real(d));

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

lambda_2 = repmat(lambda_max,3,1);

%numerical simulations

Ns = 20000; %number of realizations

theta_final = nan(Ns, N);
rr = nan(Ns, 1);

for ii_s = 1:Ns

    %random initial conditions
    rng(ii_s)
    theta0 = 2*pi*( rand(N,1) - 0.5 );

    [theta_final(ii_s,:), rr(ii_s)] = simulate_KM_final_state(a, omega, epsilon, theta0, t, dt, phi);

end

str = sprintf('files_scans/basin_size_data/N_%d_k_%d.mat', N, k);
save(str);
