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

A_wrong_indeces = [R      Q           -R*D3F2_big * 1i /phi * K_inv_tail
                   P       M           -P*D3F2_big * 1i /phi * K_inv_tail
                   Zero2   Zero1       1i / phi * K_inv_tail];
               
small_length = length(x);
n_scalar = x.size_scalar;
n_vector = x.size_vector;

big_length = size(A_wrong_indeces,1);
big_length_nodes = (big_length - n_scalar)/n_vector;
long_nodes = (big_length_nodes - 1)/2;
if mod(big_length_nodes,1)~=0 || mod(long_nodes,1)~=0
    error('Length of long vector cannot be not-integer')
end
new_nodes = long_nodes - x.nodes;
right_indices_nodes = zeros(big_length_nodes*n_vector, 1);
for j = 1: x.size_vector
    indices_j = n_scalar + (1:big_length_nodes) + (j-1)*big_length_nodes;
    first_diag = small_length + (1:new_nodes) + (2*new_nodes) * (j-1);
    second_diag = small_length + new_nodes + (1:new_nodes) + (2*new_nodes) * (j-1);
    center = n_scalar + (1:2*x.nodes+1) + (2*x.nodes+1)*(j-1);
    right_indices_nodes(indices_j) = [first_diag, center, second_diag];
end
right_indices_nodes(1:n_scalar) = 1:n_scalar;
A = A_wrong_indeces(right_indices_nodes,right_indices_nodes);


