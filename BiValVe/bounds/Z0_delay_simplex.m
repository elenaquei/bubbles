function Z0 = Z0_delay_simplex(alpha0, alpha1, alpha2, x0, x1, x2,A0_struct, A1_struct, A2_struct, DH0, DH1, DH2)

A_small0 = A0_struct.A;
A_small1 = A1_struct.A;
A_small2 = A2_struct.A;

% [A_small0, M0, P0, Q0, R0, phi0, D3F20] = A_delay_symplex(alpha0, x0);
% [A_small1, M1, P1, Q1, R1, phi1, D3F21] = A_delay_symplex(alpha1, x1);
% [A_small2, M2, P2, Q2, R2, phi2, D3F22] = A_delay_symplex(alpha2, x2);

A_small_int = interpolation(A_small0, A_small1, A_small2);

% small_der0 = derivative(alpha0, x0, 0);
% small_mat_der0 = derivative_to_matrix(small_der0);
% small_der1 = derivative(alpha1, x1, 0);
% small_mat_der1 = derivative_to_matrix(small_der1);
% small_der2 = derivative(alpha2, x2, 0);
% small_mat_der2 = derivative_to_matrix(small_der2);
% small_mat_der_int = interpolation(small_mat_der0, small_mat_der1, small_mat_der2);

small_mat_der_int = interpolation(DH0, DH1, DH2);

Z0_mat = intval(eye(size(A_small0))) - A_small_int * small_mat_der_int;


Z0 = sup(norm_Ximat(Z0_mat, x0));

end