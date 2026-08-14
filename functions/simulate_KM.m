function [theta r] = simulate_KM( w, omega, epsilon, theta0, time, dt, phi )

N = size(w,1);
theta = nan( length(time), N ); 
theta(1,:) = theta0;

r = nan(length(time),1);

r(1) = abs( (1/N) * sum( exp(1i*theta0'), 2 ) );

for ii = 2:length(time)
    
    previous_state = theta(ii-1,:);
    
    for jj = 1:N
        dth = omega(jj) + epsilon * nansum( w(jj,:) .* sin( previous_state - previous_state(jj) - phi(jj,:) ), 2 );        
        theta(ii,jj) = previous_state(jj) + (dth * dt);
    end
    
    r(ii) = abs( (1/N) * sum( exp(1i*theta(ii,:)), 2 ) );

end

% wrap theta into [-pi, pi]
theta = angle( exp( 1i*theta ) );
