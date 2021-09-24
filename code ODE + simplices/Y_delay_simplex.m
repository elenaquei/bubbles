function Y = Y_delay_simplex(alpha, x0, x1, x2)

n_nodes_long = alpha.vector_field.deg_vector*x0.nodes;

x0_long_vec = Xi_vec2vec(reshape(x0,n_nodes_long));
x1_long_vec = Xi_vec2vec(reshape(x1,n_nodes_long));
x2_long_vec = Xi_vec2vec(reshape(x2,n_nodes_long));

x_vec_int = three_intval(x0_long_vec, x1_long_vec, x2_long_vec);

Fx0 = Xi_vec2vec(apply(alpha, x0, 1));
Fx1 = Xi_vec2vec(apply(alpha, x1, 1));
Fx2 = Xi_vec2vec(apply(alpha, x2, 1));

Fx_int = three_intval(Fx0, Fx1, Fx2);

big_A0 = create_A_of_size(alpha, x0, n_nodes_long);
big_A1 = create_A_of_size(alpha, x1, n_nodes_long);
big_A2 = create_A_of_size(alpha, x2, n_nodes_long);

big_A_int = three_intval(big_A0, big_A1, big_A2);

Y_vec = x_vec_int - big_A_int * Fx_int;

Y_x_vec = vec2Xi_vec(Y_vec, x0.size_scalar, x0.size_vector, n_nodes_long);

Y = norm(Y_x_vec);

end
function val_intval = three_intval(a,b,c)

val_intval = three_intval_real(real(a), real(b), real(c)) + 1i* three_intval_real(imag(a), imag(b), imag(c));
end

function val_intval = three_intval_real(a,b,c)
val_min = min(a, min(b,c));
val_max = max(a, max(b,c));
if isintval(val_min)
    val_min = inf(val_min);
end
if isintval(val_max)
    val_max = sup(val_max);
end
val_intval = infsup(val_min, val_max);
end