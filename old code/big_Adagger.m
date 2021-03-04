function A_dagger_big = big_Adagger(A_dagger,x_PDE)

l1 = x_PDE.parameters(1);
l2 = x_PDE.parameters(2);
A_dagger_op = compose_overoperator(A_dagger,x_PDE);
A_dagger_center = reshape(A_dagger_op, 2*x_PDE.node_space, 2*x_PDE.node_time);
A_dagger_lin = linear_part_A_dagger(l1,l2,2*x_PDE.node_space, 2*x_PDE.node_time);
tail_A_dagger = tail(A_dagger_lin,x_PDE.node_space, x_PDE.node_time);
A_dagger_big = overop2mat(A_dagger_center+tail_A_dagger);