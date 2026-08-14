function A = ring_graph( N, k )

A = false( N, N );
m = -k:k;m( m == 0 ) = [];

% generate ring graph
for ii = 0:(N-1)    
    for jj = m
        
        A( ii+1, mod(ii+jj,N)+1 ) = 1;
        
    end
end
