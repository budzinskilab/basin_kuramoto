function lambda = linear_stability( A, epsilon, phi )

N = size(A,1); %number of nodes in the network

%creating matrix K
K = epsilon.*exp(-1i*phi).*A;

%eigensystem of K
[v,d] = circulant_eigensystem( K ); 

ll = repmat(real(diag(d)),3,1); %reapting eigenvalues to consider the N-periodicity of the network

aux = nan(N,N);

lambda_max = nan(N,1);

for ii = N+1:2*N
    for jj = 1:N-1
        aux(ii-N,jj) = ( ll(ii-jj) + ll(ii+jj) )./2 - ll(ii);
    end
    lambda_max(ii-N) = max(aux(ii-N,:));
end

%for plotting purposes 
aux_2 = repmat(lambda_max,1,2); 
lambda = aux_2( ceil(N/2)+1: N + ceil(N/2));
