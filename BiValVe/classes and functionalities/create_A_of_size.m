function [A,M,P,Q,R,phi, D3F2] = create_A_of_size(alpha, x, final_n_nodes, Q, R, phi, D3F2, f_n_nodes)
% can be called either with
%       function A = create_A_of_size(alpha, x, final_n_nodes)
% or with 
%       function A = create_A_of_size(A_small, M, P, Q, R, phi, D3F2, final_n_nodes)
% or with
%       function A = create_A_of_size(A_struct, final_n_nodes)
% where A_struct is a strcture with properties A, M, P, Q, R, phi, D3F2
%
% returns A, M, P, Q, R, phi, D3F2
if nargin ==2
    A_struct = alpha;
    final_n_nodes = x;
    A_small = A_struct.A;
    M = A_struct.M;
    P = A_struct.P;
    Q = A_struct.Q;
    R = A_struct.R;
    phi = A_struct.phi;
    D3F2 = A_struct.D3F2;
    size_scalar = size(R,1);
    size_vector = size(D3F2,2);
    nodes = (size(M,1)/size_vector - 1)/2;
    small_length = size(A_small,1);
    if final_n_nodes < nodes
        error('Cannot have a matrix A of size smaller than the initial vector')
    end
elseif nargin ==3
    if final_n_nodes < x.nodes
        error('Cannot have a matrix A of size smaller than the initial vector')
    end
    [A_small, M, P, Q, R, phi, D3F2] = A_delay_symplex(alpha, x);
    nodes = x.nodes;
    size_vector = x.size_vector;
    small_length = length(x);
elseif nargin == 8
    A_small = alpha;
    M = x;
    P = final_n_nodes;
    final_n_nodes = f_n_nodes;
    small_length = size(A_small,1);
    size_scalar = size(R,1);
    size_vector = size(D3F2,2);
    nodes = ((small_length-size_scalar) / size_vector-1)/2;
end

if final_n_nodes == nodes
    A = A_small;
    return
end

K_tail = [-final_n_nodes:-nodes-1, nodes+1:final_n_nodes];

K_inv_tail = diag(1./repmat(K_tail,1,size_vector));

Zero1 = zeros(size(K_inv_tail,1), size(Q,2));
Zero2 = zeros(size(K_inv_tail,1), size(R,2));

D3F2_big = kron(D3F2, ones(1,length(K_tail)));

A_wrong_indeces = [R      Q           -R*D3F2_big * 1i /phi * K_inv_tail
                   P       M           -P*D3F2_big * 1i /phi * K_inv_tail
                   Zero2   Zero1       1i / phi * K_inv_tail];

big_length = size(A_wrong_indeces,1);
big_length_nodes = (big_length - size_scalar)/size_vector;
long_nodes = (big_length_nodes - 1)/2;
if mod(big_length_nodes,1)~=0 || mod(long_nodes,1)~=0
    error('Length of long vector cannot be not-integer')
end
new_nodes = long_nodes - nodes;
right_indices_nodes = zeros(big_length_nodes*size_vector, 1);
for j = 1: size_vector
    indices_j = size_scalar + (1:big_length_nodes) + (j-1)*big_length_nodes;
    first_diag = small_length + (1:new_nodes) + (2*new_nodes) * (j-1);
    second_diag = small_length + new_nodes + (1:new_nodes) + (2*new_nodes) * (j-1);
    center = size_scalar + (1:2*nodes+1) + (2*nodes+1)*(j-1);
    right_indices_nodes(indices_j) = [first_diag, center, second_diag];
end
right_indices_nodes(1:size_scalar) = 1:size_scalar;
A = A_wrong_indeces(right_indices_nodes,right_indices_nodes);

if nargout <2
    return
end

index_scalar = 1:size_scalar;
index_vector = size_scalar+1:size(A,1);

M = A(index_vector,index_vector);
P = A(index_vector,index_scalar);
Q = A(index_scalar, index_vector);
R = A(index_scalar, index_scalar);



