
function phi = heterogeneous_phase_lag( a, max_phi, k )

N = size( a, 1 ); phi = zeros( N );
for ii = 1:N
    for jj = 1:N
        if ( a(ii,jj) > 0 )
            d = abs( ii - jj ); d = min( N - d, d ); phi(ii,jj) = (d/k)*(max_phi);
        else
            phi(ii,jj) = 0;
        end
    end
end