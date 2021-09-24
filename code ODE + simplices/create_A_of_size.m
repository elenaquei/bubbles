function A = create_A_of_size(alpha, x, final_n_nodes)

if final_n_nodes < x.nodes
    error('Cannot have a matrix A of size smaller than the initial vector')
end

[A_small, M, P, Q, R, phi, D3F2] = A_delay_symplex(alpha, x);

if final_n_nodes == x.nodes
    A = A_small;
    return
end

K_tail = [-final_n_nodes:-x.nodes-1, x.nodes+1:final_n_nodes];

K_inv_tail = diag(1./repmat(K_tail,1,x.size_vector));

Zero1 = zeros(size(K_inv_tail,1), size(Q,2));
Zero2 = zeros(size(K_inv_tail,1), size(R,2));

D3F2_big = kron(D3F2, ones(1,length(K_tail)));

A = [M     P      -P*D3F2_big * 1i /phi * K_inv_tail
    Q      R      -R*D3F2_big * 1i /phi * K_inv_tail
    Zero1  Zero2  1i / phi * K_inv_tail];

