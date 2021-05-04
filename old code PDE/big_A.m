function A_big = big_A(A,x_PDE,new_dim_PDE)

if nargin==2
    new_dim_PDE = reshape(x_PDE,2*x_PDE.node_space, 2*x_PDE.node_time);
end
if size(A) == length(new_dim_PDE)
    A_big = A;
    return
end


l1 = x_PDE.parameters(1);
l2 = x_PDE.parameters(2);
A_op = compose_overoperator(A,x_PDE);
A_center = reshape(A_op, new_dim_PDE.node_space, new_dim_PDE.node_time);
A_lin = linear_part_A(l1,l2,new_dim_PDE.node_space, new_dim_PDE.node_time);
tail_A = tail(A_lin,x_PDE.node_space, x_PDE.node_time);
A_big = overop2mat(A_center+tail_A);