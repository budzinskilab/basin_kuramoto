%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                            %%
%% Linear stability and basin size estimatation in Kuramoto networks          %%
%% Roberto Budzinski & Kalel Rossi                                            %%
%% August 2026                                                                %%
%%                                                                            %%
%% Figure 1 — plotting estimatation and simulations of basin size of q-states %%
%%                                                                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% setup

clearvars; clc;
addpath( './functions' );

%% examples — basin size estimation

%example 1

ii_jj = 6;

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

cc = cc ./ max(cc);
cc = cc ./ nansum(cc);

%estimating basin size
gamma = diag( real(d) );
for ii = 1:N
    if ( lambda_max(ii) >= 0 )
        gamma(ii) = 0;
    end
end

gamma = power( gamma, N );
gamma = gamma ./ nansum(gamma);

gamma( gamma == 0 ) = nan;
cc( cc == 0 ) = nan;

aux = repmat(gamma,1,2);
gamma = aux( ceil(N/2)+1: N + ceil(N/2));

aux = repmat(cc,1,2);
cc = aux( ceil(N/2)+1: N + ceil(N/2));

fg1 = figure; hold on
stem(cc, 'filled', 'MarkerSize', 10, 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5]);
stem(gamma, 'filled', 'Color', [0 0 0],  'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
ylabel('relative basin size')
xlabel('q-state');
lb = legend('numerical', 'prediction');
lb.EdgeColor = 'white';
lb.Location = 'north';
ylim([-0.05 0.55]);
yticks([0 0.25 0.5]);
xlim([0 52])
xticks([6 16 26 36 46])
xticklabels({'-20', '-10', '0', '10', '20'});
set(gcf,'Position', [1281         747         334         221]);
set( gca, 'fontname', 'arial', 'fontsize', 13, 'linewidth', 2 );
box off

% example 2

ii_jj = 15;

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

cc = cc ./ max(cc);
cc = cc ./ nansum(cc);

%estimating basin size
gamma = diag( real(d) );
for ii = 1:N
    if ( lambda_max(ii) >= 0 )
        gamma(ii) = 0;
    end
end

gamma = power( gamma, N );
gamma = gamma ./ nansum(gamma);

gamma( gamma == 0 ) = nan;
cc( cc == 0 ) = nan;

aux = repmat(gamma,1,2);
gamma = aux( ceil(N/2)+1: N + ceil(N/2));

aux = repmat(cc,1,2);
cc = aux( ceil(N/2)+1: N + ceil(N/2));

fg2 = figure; hold on
stem(cc, 'filled', 'MarkerSize', 10, 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5]);
stem(gamma, 'filled', 'Color', [0 0 0],  'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
ylabel('relative basin size')
xlabel('q-state');
lb = legend('numerical', 'prediction');
lb.EdgeColor = 'white';
lb.Location = 'north';
ylim([-0.05 0.55]);
yticks([0 0.25 0.5]);
xlim([0 52])
xticks([6 16 26 36 46])
xticklabels({'-20', '-10', '0', '10', '20'});
set(gcf,'Position', [ 1000         748         334         224]);
set( gca, 'fontname', 'arial', 'fontsize', 13, 'linewidth', 2 );
box off

%example 3

ii_jj = 19;

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

cc = cc ./ max(cc);
cc = cc ./ nansum(cc);

%estimating basin size
gamma = diag( real(d) );
for ii = 1:N
    if ( lambda_max(ii) >= 0 )
        gamma(ii) = 0;
    end
end

gamma = power( gamma, N );
gamma = gamma ./ nansum(gamma);

gamma( gamma == 0 ) = nan;
cc( cc == 0 ) = nan;

aux = repmat(gamma,1,2);
gamma = aux( ceil(N/2)+1: N + ceil(N/2));

aux = repmat(cc,1,2);
cc = aux( ceil(N/2)+1: N + ceil(N/2));


fg3 = figure; hold on
stem(cc, 'filled', 'MarkerSize', 10, 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5]);
stem(gamma, 'filled', 'Color', [0 0 0],  'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0]);
ylabel('relative basin size')
xlabel('q-state');
lb = legend('numerical', 'prediction');
lb.EdgeColor = 'white';
lb.Location = 'north';
ylim([-0.05 0.55]);
yticks([0 0.25 0.5]);
xlim([0 52])
xticks([6 16 26 36 46])
xticklabels({'-20', '-10', '0', '10', '20'});
set(gcf,'Position', [712   743   334   223]);
set( gca, 'fontname', 'arial', 'fontsize', 13, 'linewidth', 2 );
box off

%% error analysis — basin size as a function of k

N = 51;

ee = nan( N, floor(N/2) );
M = nan( floor(N/2), 1 );
error_basin = nan( floor(N/2), 1 );

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

    cc = cc./max(cc);
    cc = cc ./ nansum(cc);

    %estimating basin size
    gamma = diag( real(d) );
    for ii = 1:N
        if ( lambda_max(ii) >= 0 )
            gamma(ii) = 0;
        end
    end

    gamma = power( gamma, N );
    gamma = gamma ./ nansum(gamma);

    gamma( gamma == 0 ) = nan;
    cc( cc == 0 ) = nan;

    aux = 0;
    for ii = 1:N
        if ( lambda_max(ii) < 0 )
            ee(ii,ii_jj) = abs( cc(ii) - gamma(ii) );
            aux = aux +1;
        end
    end

    M(ii_jj) = aux;

    error_basin(ii_jj) = nansum( ee(:,ii_jj) );
    error_basin(ii_jj) = error_basin(ii_jj) ./ M(ii_jj);

end

% Plot the error basin as a function of k
fg4 = figure;
bar(1:floor(N/2), error_basin);
ylim( [-0.01 0.3] );
xlabel('k');
ylabel('Mean error — basin size estimation');
set(gcf,'Position', [734   302   834   373]);
set( gca, 'fontname', 'arial', 'fontsize', 14, 'linewidth', 2 );
box off