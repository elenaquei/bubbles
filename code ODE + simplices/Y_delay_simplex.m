function Y = Y_delay_simplex(alpha0, alpha1, alpha2, x0, x1, x2)
alpha = alpha0;
n_nodes_long = alpha.vector_field.deg_vector*x0.nodes;

x0_long_vec = Xi_vec2vec(reshape(x0,n_nodes_long));
x1_long_vec = Xi_vec2vec(reshape(x1,n_nodes_long));
x2_long_vec = Xi_vec2vec(reshape(x2,n_nodes_long));

x_vec_int = interpolation(x0_long_vec, x1_long_vec, x2_long_vec);

x_int = vec2Xi_vec(x_vec_int, x0.size_scalar, x0.size_vector, n_nodes_long);

alpha  = interpolation(alpha0, alpha1, alpha2);
Fx_int = Xi_vec2vec(apply(alpha, x_int));

big_A0 = create_A_of_size(alpha0, x0, n_nodes_long);
big_A1 = create_A_of_size(alpha1, x1, n_nodes_long);
big_A2 = create_A_of_size(alpha2, x2, n_nodes_long);

big_A_int = interpolation(big_A0, big_A1, big_A2);

Y_vec = big_A_int * Fx_int;

Y_x_vec = vec2Xi_vec(Y_vec, x0.size_scalar, x0.size_vector, n_nodes_long);

Y = norm(Y_x_vec);

end