%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                           %%
%% Linear stability and basin size estimatation in Kuramoto networks         %%
%% Roberto Budzinski & Kalel Rossi                                           %%
%% August 2026                                                               %%
%%                                                                           %%
%% Figure 1 — calculation and plot of basin geometry                         %%
%%                                                                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% setup
clearvars; clc;
addpath( './functions' );

%% parameters

dt = 0.001; T = 10.0; t = 0:dt:T;   %time
N = 51;                             %number of oscillators
epsilon = 10;                       %coupling strength
f_mu = 0;                           %(Hz) natural frequency
omega = ( f_mu*ones(N,1) )*2*pi;    %natural frequency

Ns = 1000;

%% example 1: k = 15

ii_k = 15;

% adjacency matrix
k = ii_k;
A = double( ring_graph( N, k ) );

% phase-lags
phi = heterogeneous_phase_lag( A, pi, k );

%composite matrix
K = epsilon.*exp(-1i*phi).*A;

%eigenspectrum of K
[v,~] = circulant_eigensystem( K ); %obtaing eigenvectors of K

%initial condition — synchronized state
theta0 = zeros(N,1);

%perturbations
alpha1 = linspace(-pi, pi, Ns);
alpha2 = linspace(-pi, pi, Ns);

%similarity measure
s = nan(Ns);

for ii_s = 1:Ns
    for jj_s = 1:Ns

        theta0_pert = theta0;

        theta0_pert(1) = theta0(1) + alpha1(ii_s);
        theta0_pert(2) = theta0(2) + alpha2(jj_s);

        %numerical simulation
        theta_final = simulate_KM_final_state(A, omega, epsilon, theta0_pert, t, dt, phi);

        for ii = 1:N
            s_aux = order_parameter( theta_final - angle(v(:,ii) )', N );
            if( s_aux > 0.99)
                if( ii <= floor(N/2))
                    ii_aux = ii - 1;
                else
                    ii_aux = ii - (N + 1);
                end
                s(ii_s, jj_s) = ii_aux;
            end
        end
    end
end

%plot

rng(31)
rcm = rand(5,3);
cm = colormap(rcm);

fg1 = figure;
imagesc(alpha1, alpha2, s);
xlabel('\alpha_1');
ylabel('\alpha_2');
xticks([-pi, 0, pi]);
yticks([-pi, 0, pi]);
xticklabels({'-\pi', '0', '\pi'});
yticklabels({'-\pi', '0', '\pi'});
colormap(cm);
set( gca, 'fontname', 'arial', 'fontsize', 15, 'linewidth', 2 );
set(gcf,'Position', [110   562   251   216]);

%% example 2: k = 19

ii_k = 19;

% adjacency matrix
k = ii_k;
A = double( ring_graph( N, k ) );

% phase-lags
phi = heterogeneous_phase_lag( A, pi, k );

%composite matrix
K = epsilon.*exp(-1i*phi).*A;

%eigenspectrum of K
[v,~] = circulant_eigensystem( K ); %obtaing eigenvectors of K

%initial condition — synchronized state
theta0 = zeros(N,1);

%perturbations
alpha1 = linspace(-pi, pi, Ns);
alpha2 = linspace(-pi, pi, Ns);

%similarity measure
s = nan(Ns);

for ii_s = 1:Ns
    for jj_s = 1:Ns

        theta0_pert = theta0;

        theta0_pert(1) = theta0(1) + alpha1(ii_s);
        theta0_pert(2) = theta0(2) + alpha2(jj_s);

        %numerical simulation
        theta_final = simulate_KM_final_state(A, omega, epsilon, theta0_pert, t, dt, phi);

        for ii = 1:N
            s_aux = order_parameter( theta_final - angle(v(:,ii) )', N );
            if( s_aux > 0.99)
                if( ii <= floor(N/2))
                    ii_aux = ii - 1;
                else
                    ii_aux = ii - (N + 1);
                end
                s(ii_s, jj_s) = ii_aux;
            end
        end
    end
end

%plot

rng(31)
rcm = rand(5,3);
cm = colormap(rcm);

fg2 = figure;
imagesc(alpha1, alpha2, s);
xlabel('\alpha_1');
ylabel('\alpha_2');
xticks([-pi, 0, pi]);
yticks([-pi, 0, pi]);
xticklabels({'-\pi', '0', '\pi'});
yticklabels({'-\pi', '0', '\pi'});
colormap(cm);
colorbar;
set( gca, 'fontname', 'arial', 'fontsize', 15, 'linewidth', 2 );
set(gcf,'Position', [110   562   251   216]);


