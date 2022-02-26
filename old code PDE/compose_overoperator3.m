function A_op= compose_overoperator3(A, size_real, node_space, node_time)
% function A_op= compose_overoperator(A, size_real, node_space, node_time)
%
% INPUT 
% A             big matrix 
% size_real     number of complex variables FOR ONE VECTOR
% node_space    number of modes used for space discretization FOR ONE
% VECTOR
% node_time     number of nodes used for time discretization
% OUTPUT
% A_op          overoperator correspongin to the matrix A
% 
% or
%
% INPUT
% A             big matrix
% X_PDE_vec     PDE_vector of the wanted dimension

if nargin == 2
    X_PDE = size_real;
    size_real = X_PDE.size_real;
    node_space = X_PDE.node_space;
    node_time = X_PDE.node_time;
end
size_real = size_real/3;

real_indices = 1:size_real;
l1_indices = real_indices(end) + (1:2*node_space+1);
L1_indices = l1_indices(end) + (1:(2*node_space+1)*(2*node_time+1));

index_1 = 1:L1_indices(end);
index_2 = L1_indices(end) + index_1;
index_3 = L1_indices(end)*2 +index_1;

A11 = A(index_1,index_1);
A11 = compose_overoperator(A11, size_real, node_space, node_time);
A12 = A(index_1, index_2);
A12 = compose_overoperator(A12, size_real, node_space, node_time);
A13 = A(index_1, index_3);
A13 = compose_overoperator(A13, size_real, node_space, node_time);

A21 = A(index_2,index_1);
A21 = compose_overoperator(A21, size_real, node_space, node_time);
A22 = A(index_2, index_2);
A22 = compose_overoperator(A22, size_real, node_space, node_time);
A23 = A(index_2, index_3);
A23 = compose_overoperator(A23, size_real, node_space, node_time);

A31 = A(index_3,index_1);
A31 = compose_overoperator(A31, size_real, node_space, node_time);
A32 = A(index_3, index_2);
A32 = compose_overoperator(A32, size_real, node_space, node_time);
A33 = A(index_3, index_3);
A33 = compose_overoperator(A33, size_real, node_space, node_time);

A_op = overoperator3(A11, A12, A13, A21, A22, A23, A31, A32, A33);

