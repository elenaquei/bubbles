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

norm_P = norm_ell1_to_C(P_int);
norm_R = norm_C_to_C(R_int);

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
% missing D_3F^2*i*psi!!

% STAR2
Pi_infty_indices = n_nodes*alpha.deg_vector + 1 + ...
    [-n_nodes*alpha.deg_vector:-n_nodes-1, n_nodes+1:n_nodes*alpha.deg_vector];
derivative_parameters_F = DF_derivative.derivative_Flambda;
derivative_parameters_F_tail = derivative_parameters_F.*(Pi_infty_indices~=0);
STAR2 = Xi_vector([],derivative_parameters_F_tail);
norm_STAR2 = norm(STAR2);
norm_Kinv_STAR2 = sum(1./Pi_infty_indices) * norm_STAR2;

norm_A_STAR2 = [norm_P * norm_Kinv_STAR2
    norm_R * norm_Kinv_STAR2
    abs(phi_int)^-1 * norm_Kinv_STAR2];

% STAR3
norm_STAR3 = zeros(alpha.vector_field.n_equations);
for i =1:N
    for j = 1:N
        norm_STAR3(i,j) = operator_tail_norm(A.derivative_Fx_toeplix(j,i), ...
            nodes, size_vector, alpha.deg_vector);
        delay_part = 0;
        for term = 1:size(A.derivative_Fx_delay{i,j})
            if ~isempty(A.derivative_Fx_delay{i,j}{term})
                conv = A.derivative_Fx_delay{i,j}{term}.convolution;
                norm_conv = operator_tail_norm(conv, nodes, size_vector, ...
                    alpha.deg_vector);
                delay_part = delay_part +  norm_conv;
            end
        end
        norm_STAR3(i,j) = norm_STAR3(i,j) + delay_part;
    end
end
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

function n = norm_ell1_to_C(P, nodes, size_scalar, size_vector)
n = zeros(size_vector,1);
K = (-nodes:nodes).';
K_big = repmat(K, 1, size_scalar);
for i = 1:size_vector
    indeces = (1:2*nodes+1) + (i-1) * (2*nodes+1);
    n(i) = sum(abs(P(indeces, :)).* nu.^(-abs(K_big)),2);
end
end

function n = norm_C_to_C(R)
n = zeros(size(R,1));
n = sum(abs(R),2);
end

function operator_tail_norm(x, nodes, size_vector, degree)
norm_op = size(size_vector);
for i = 1:size_vector
    index_i = 1:2*nodes+1 + (i-1) * (2*nodes+1);
    for j = 1:size_vector
        index_j = 1:2*nodes+1 + (j-1) * (2*nodes+1);
        norm_op(i,j) = sub_operator_tail_norm(x(index_i, index_j), ...
            nodes, degree);
    end
end

end

function sub_operator_tail_norm(x, nodes, degree)
K_all = -nodes*degree : nodes*degree;
K_center = K_all * (abs(K_all)<=nodes);
K_tail = K_all * (abs(K_all)>nodes);

% max |m|>nodes sum |n|<nodes |x_{n-m}| * nu ^(|n|-|m|)
% if n-m>deg*nodes, x_{n-m} = 0

end