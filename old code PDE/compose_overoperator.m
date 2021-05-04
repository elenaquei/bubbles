function A_op= compose_overoperator(A, size_real, node_space, node_time)
% function A_op= compose_overoperator(A, size_real, node_space, node_time)
%
% INPUT 
% A             big matrix 
% size_real     number of complex variables
% node_space    number of modes used for space discretization
% node_time     number of nodes used for time discretization
% OR
% A             big matrix 
% X             PDE_vector
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

real_indices = 1:size_real;
l1_indices = real_indices(end) + (1:2*node_space+1);
L1_indices = l1_indices(end) + (1:(2*node_space+1)*(2*node_time+1));

A11 = A(real_indices,real_indices);
A12 = A(real_indices, l1_indices);
A13 = A(real_indices, L1_indices);

A21 = A(l1_indices,real_indices);
A22 = A(l1_indices, l1_indices);
A23 = A(l1_indices, L1_indices);

A31 = A(L1_indices,real_indices);
A32 = A(L1_indices, l1_indices);
A33 = A(L1_indices, L1_indices);

A_op = overoperator(A11, A12, A13, A21, A22, A23, A31, A32, A33);