function [Z2vector,Z2_s]=Z2_delay_simplex(x0,x1,x2,alpha, RAD_MAX)
global nu

M = x0.n_nodes;

[A_small0, M0, P0, Q0, R0, phi0, D3F20] = A_delay_symplex(alpha, x0);
[A_small1, M1, P1, Q1, R1, phi1, D3F21] = A_delay_symplex(alpha, x1);
[A_small2, M2, P2, Q2, R2, phi2, D3F22] = A_delay_symplex(alpha, x2);

A_small_int = three_intval(A_small0, A_small1, A_small2);

M_int = three_intval(M0, M1, M2);
P_int = three_intval(P0, P1, P2);
Q_int = three_intval(Q0, Q1, Q2);
R_int = three_intval(R0, R1, R2);
phi_int = three_intval(phi0, phi1, phi2);
D3F2_int = three_intval(D3F20, D3F21, D3F22);

block_norm_A = [1,2,3]; %TODO

upper_bound_psi_DF = [1,2,3]; %TODO

Z22 = block_norm_A * upper_bound_psi_DF;

if nu^(2*M+2) < exp(1)
    Lemma_bound = nu^(2*M+2)/ (exp(1) * log(nu^(2*M+2)));
else
    Lemma_bound = 1;
end

block_norm_AM = [1,2,3];
upper_bound_DDF_without_dpsi = [1,2,3];
Z21 = block_norm_AM * upper_bound_DDF_without_dpsi;

Z2 = Z22 + Z21;
end

function block_norm(M, P, Q, R, V)

end

function upper_bound_psi_DF()
end

function upper_bound_DDF_without_psi()
end