function Z0 = Z0_delay_simplex(alpha0, alpha1, alpha2, x0, x1, x2)
global use_intlab
use_intlab = 0;
[A_small0, M0, P0, Q0, R0, phi0, D3F20] = A_delay_symplex(alpha0, x0);
[A_small1, M1, P1, Q1, R1, phi1, D3F21] = A_delay_symplex(alpha1, x1);
[A_small2, M2, P2, Q2, R2, phi2, D3F22] = A_delay_symplex(alpha2, x2);

A_small_int = three_intval(A_small0, A_small1, A_small2);

% M_int = three_intval(M0, M1, M2);
% P_int = three_intval(P0, P1, P2);
% Q_int = three_intval(Q0, Q1, Q2);
% R_int = three_intval(R0, R1, R2);
% phi_int = three_intval(phi0, phi1, phi2);
% D3F2_int = three_intval(D3F20, D3F21, D3F22);


use_intlab = 1;
small_der0 = derivative(alpha0, x0, 0);
small_mat_der0 = derivative_to_matrix(small_der0);
small_der1 = derivative(alpha1, x1, 0);
small_mat_der1 = derivative_to_matrix(small_der1);
small_der2 = derivative(alpha2, x2, 0);
small_mat_der2 = derivative_to_matrix(small_der2);

small_mat_der_int = three_intval(small_mat_der0, small_mat_der1, small_mat_der2);

Z0_mat = intval(eye(size(A_small0))) - A_small_int * small_mat_der_int;


Z0 = norm_Ximat(Z0_mat, x0);

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