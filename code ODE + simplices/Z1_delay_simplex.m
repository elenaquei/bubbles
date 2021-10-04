function Z1 = Z1_delay_simplex(alpha, x0, x1, x2)

x0_long_vec = Xi_vec2vec(reshape(x0,n_nodes_long));
x1_long_vec = Xi_vec2vec(reshape(x1,n_nodes_long));
x2_long_vec = Xi_vec2vec(reshape(x2,n_nodes_long));

x_vec_int = three_intval(x0_long_vec, x1_long_vec, x2_long_vec);

x_int = vec2Xi_vec(x_vec_int, x0.size_scalar, x0.size_vector, n_nodes_long);

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

DF_derivative = derivative(alpha, x);

norm_P = norm_ell1_to_C(P);
norm_R = norm_C_to_C(R);

% first term: A_s * Pi_M ( D_z Pi_infty F)
n_nodes = x0.nodes;

% norm of STAR1
norm_STAR1 = zeros(alpha.vector_field.n_equations);
N = alpha.vector_field.n_equations;
for i =1:N
    for j = 1:N
        norm_STAR1(i,j) = norm(Xi_vector([],A.derivative_Fx_toeplix(j,i)));
        delay_part = 0;
        for term = 1:size(A.derivative_Fx_delay{i,j})
            if ~isempty(A.derivative_Fx_delay{i,j}{term})
                conv = A.derivative_Fx_delay{i,j}{term}.convolution;
                norm_conv = norm(Xi_vector([], conv));
                delay_part = delay_part +  norm_conv;
            end
        end
        norm_STAR1(i,j) = norm_STAR1(i,j) + delay_part;
    end
end
norm_Kinv_STAR1 = 1/(n_nodes+1) * norm_STAR1;

norm_A_STAR1 = [norm_P * norm_Kinv_STAR1
    norm_R * norm_Kinv_STAR1
    abs(phi_int)^-1 * norm_Kinv_STAR1];

% STAR2
Pi_infty_indices = n_nodes_long
STAR2 = derivative_Flambda;
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

function n = norm_ell1_to_C(P)
n = zeros(size(P,1));
nodes = (zeros(size(P,2))-1)/2;
K = -nodes:nodes;
K_big = repmat(K, n, 1);
n = sum(abs(P).* nu.^(-abs(K_big)),2);
end

function n = norm_C_to_C(R)
n = zeros(size(R,1));
n = sum(abs(R),2);
end