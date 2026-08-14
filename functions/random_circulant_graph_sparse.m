function [W, gen_vec] = random_circulant_graph_sparse( N, M, kk )

rng(0)
gen_vec = rand(M,1);

W_vec = zeros(N,1);

rng(kk);
idx = randperm(N);
idx = idx(1:M);

for jj = 1:M
    W_vec(idx(jj)) = gen_vec(jj);
end

W = zeros(N);

for i = 1:N
    W(i,:) = W_vec([N-i+1:N 1:N-i]);
end
