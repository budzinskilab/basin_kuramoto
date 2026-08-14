function [theta r] = simulate_KM_final_state( w, omega, epsilon, theta0, time, dt, phi )

N = size(w,1);
theta = theta0';

for ii = 2:length(time)
    
    previous_state = theta;
    
    for jj = 1:N
        dth = omega(jj) + epsilon * nansum( w(jj,:) .* sin( previous_state - previous_state(jj) - phi(jj,:) ), 2 );
        theta(jj) = previous_state(jj) + (dth * dt);
    end

    theta = angle(exp(1i*theta));
    r = abs( (1/N) * sum( exp(1i*theta), 2 ) );
       
end
