%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                   %%
%% Linear stability and basin size estimatation in Kuramoto networks %%
%% Roberto Budzinski & Kalel Rossi                                   %%
%% August 2026                                                       %%
%%                                                                   %%
%% Figure 2 — ratio of q-states over all realizations                %%
%%                                                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% setup

clearvars; clc;
addpath( './functions' );

%% parameters

N = 51;        %number of oscillators

%% basin sizes for all states
q_states = nan( floor(N/2),1 );

for ii_jj = 1:floor(N/2)

    str = sprintf('files_scans/basin_size_data/N_51_k_%d.mat', ii_jj);
    load(str);

    s  = nan(Ns,N);
    cc = zeros(N,1);

    for jj = 1:N
        for ii = 1:Ns
            s(ii,jj) = order_parameter( theta_final(ii,:) - angle(v(:,jj) )', N );
        end
        cc(jj) = nansum( s(:,jj) > 0.99 );
    end
    cc = cc./Ns;
    q_states(ii_jj) = nansum(cc);

end

fg1 = figure;
bar(q_states);
xlim([0.25 25.75]);
ylim([0 1.05]);
xticks([1 5 10 15 20 25]);
xlabel('k');
ylabel('ratio of q-states');
set(gcf,'Position', [ 3   782   802   247]);
set( gca, 'fontname', 'arial', 'fontsize', 17, 'linewidth', 2 );
box off
