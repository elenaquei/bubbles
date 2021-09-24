function [A_small, M, P, Q, R, phi, D3F2] = A_delay_symplex(alpha, x)

small_der = derivative(alpha, x, 0);
small_mat_der = derivative_to_matrix(small_der);

inverse_mat_der = inv(small_mat_der);

index_scalar = 1:x.size_scalar;
index_vector = x.size_scalar+1:length(x);

M = inverse_mat_der(index_vector,index_vector);
P = inverse_mat_der(index_vector,index_scalar);
Q = inverse_mat_der(index_scalar, index_vector);
R = inverse_mat_der(index_scalar, index_scalar);

phi = x.scalar(1);

% attention, attention, attention: D3F2, as defined, is a hell of a mess!
D3F2 = small_der.derivative_G_hat;
D3F2(1:alpha.scalar_equations.number_equations_pol) = 0;

A_small = inverse_mat_der;
