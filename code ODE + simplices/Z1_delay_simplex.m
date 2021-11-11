function Z1 = Z1_delay_simplex(alpha0, alpha1, alpha2, x0, x1, x2)
global use_intlab
use_intlab = 1;

degree = alpha0.vector_field.deg_vector;
n_nodes_long = x0.nodes*degree;

x_int = interpolation(x0, x1, x2);

[A_small0, M0, P0, Q0, R0, phi0, D3F20] = A_delay_symplex(alpha0, x0);
[A_small1, M1, P1, Q1, R1, phi1, D3F21] = A_delay_symplex(alpha1, x1);
[A_small2, M2, P2, Q2, R2, phi2, D3F22] = A_delay_symplex(alpha2, x2);

A_small_int = interpolation(A_small0, A_small1, A_small2);

M_int = interpolation(M0, M1, M2);
P_int = interpolation(P0, P1, P2);
Q_int = interpolation(Q0, Q1, Q2);
R_int = interpolation(R0, R1, R2);
phi_int = interpolation(phi0, phi1, phi2);
D3F2_int = interpolation(D3F20, D3F21, D3F22);

bool_long = 1; 
% we already have extended x to be long enough, we don't need additional
% terms
alpha_int = interpolation(alpha0, alpha1, alpha2);
DF = derivative(alpha_int, x_int, bool_long);

norm_D3F2 = intval(sup(abs(D3F2_int)));
norm_P = intval(norm_C_to_ell1(P_int, x0.nodes, x0.size_scalar, x0.size_vector));
norm_R = abs(R_int);
norm_M = intval(norm_ell1_to_ell1(M_int, x0.nodes, x0.size_vector));
norm_Q = intval(norm_ell1_to_C(Q_int, x0.nodes, x0.size_scalar, x0.size_vector));
% first term: A_s * Pi_M ( D_z Pi_infty F)
n_nodes = x0.nodes;

% norm of STAR1
norm_STAR1 = intval(zeros(alpha0.vector_field.n_equations));
N = alpha0.vector_field.n_equations;
for i =1:N
    for j = 1:N
        norm_STAR1(i,j) = intval(norm(Xi_vector([],DF.derivative_Fx_toeplix(j,i,:))));
        delay_part = 0;
        for term = 1:size(DF.derivative_Fx_delay{i,j})
            if ~isempty(DF.derivative_Fx_delay{i,j}{term})
                conv = DF.derivative_Fx_delay{i,j}{term}.convolution;
                norm_conv = intval(norm(Xi_vector([], conv.')));
                delay_part = delay_part +  norm_conv;
            end
        end
        norm_STAR1(i,j) = norm_STAR1(i,j) + delay_part;
    end
end
norm_Kinv_STAR1 = 1/intval(n_nodes+1) * norm_STAR1;

norm_A_STAR1 = [norm_P * norm_D3F2 * abs(phi_int)^-1 * norm_Kinv_STAR1
    norm_R * norm_D3F2 * abs(phi_int)^-1 * norm_Kinv_STAR1
    abs(phi_int)^-1 * norm_Kinv_STAR1];

% STAR2
Pi_infty_indices = ...
    [-n_nodes*degree:-n_nodes-1, 0*(-n_nodes:n_nodes), n_nodes+1:n_nodes*degree];
center_indices = n_nodes*degree-1:n_nodes*degree + 2*n_nodes-1;
derivative_parameters_F_tail = DF.derivative_Flambda;
norm_STAR2 = intval(zeros(x0.size_vector, x0.size_scalar));
for i = 1:x0.size_scalar
    for j = 1:x0.size_vector
        derivative_parameters_F_tail(i,j, center_indices) = 0;
        norm_STAR2(j,i) = intval(norm(Xi_vector([],derivative_parameters_F_tail(i,j,:))));
    end
end
norm_Kinv_STAR2 = sum(abs(1./intval(Pi_infty_indices(Pi_infty_indices~=0)))) * norm_STAR2;


norm_A_STAR2 = [norm_P * norm_D3F2 * abs(phi_int)^-1 * norm_Kinv_STAR2
    norm_R.' * norm_D3F2 * abs(phi_int)^-1 * norm_Kinv_STAR2
    abs(phi_int)^-1 * norm_Kinv_STAR2];

% STAR3
norm_STAR3 = intval(zeros(alpha0.vector_field.n_equations));
for i =1:N
    for j = 1:N
        norm_STAR3(i,j) = intval(operator_tail_norm(DF.derivative_Fx_toeplix(j,i,:), ...
            n_nodes, degree));
        delay_part = 0;
        for term = 1:size(DF.derivative_Fx_delay{i,j})
            if ~isempty(DF.derivative_Fx_delay{i,j}{term})
                conv = DF.derivative_Fx_delay{i,j}{term}.convolution;
                norm_conv = operator_tail_norm(conv, n_nodes, ...
                    degree);
                delay_part = delay_part +  norm_conv;
            end
        end
        norm_STAR3(i,j) = norm_STAR3(i,j) + delay_part;
    end
end
% I think??
norm_Kinv_STAR3 = 1/intval(n_nodes+1) * norm_STAR3;

norm_A_STAR3 = [norm_M * norm_Kinv_STAR3
    norm_Q * norm_Kinv_STAR3
    zeros(x0.size_vector,x0.size_vector) * norm_Kinv_STAR3];
% STAR4

warning('Ordering of rows need to be consistent!')
Z1 = sum([norm_A_STAR1, norm_A_STAR2, norm_A_STAR3 + norm_A_STAR1],2);
end

function n = norm_C_to_ell1(P, nodes, size_scalar, size_vector)
global use_intval
global nu
n = zeros(size_vector,size_scalar);
if isintval(P) || use_intval
    n = intval(n);
end
K = (-nodes:nodes).';
K_big = repmat(K, 1, size_scalar);
for i = 1:size_vector
    indeces = (1:2*nodes+1) + (i-1) * (2*nodes+1);
    n(i,:) = sum(abs(P(indeces, :)).* nu.^(abs(K_big)),1);
end
end

function n = norm_ell1_to_C(Q, nodes, size_scalar, size_vector)
global use_intval
global nu
n = zeros(size_scalar, size_vector);
if isintval(Q) || use_intval
    n = intval(n);
end
K = (-nodes:nodes);
K_big = repmat(K, size_scalar, 1);
for i = 1:size_vector
    indeces = (1:2*nodes+1) + (i-1) * (2*nodes+1);
     n(:,i) = max(abs(Q(:, indeces)).* nu.^(-abs(K_big)), [], 2);
end
end

function n = norm_ell1_to_ell1(M, nodes, size_vector)
global use_intval
global nu
n = zeros(size_vector, size_vector);
if isintval(M) || use_intval
    n = intval(n);
end
K = (-nodes:nodes);
for i = 1:size_vector
    indeces_i = (1:2*nodes+1) + (i-1) * (2*nodes+1);
    for j = 1:size_vector
        indeces_j = (1:2*nodes+1) + (j-1) * (2*nodes+1);
        n(i,j) = max(nu.^(abs(K))*abs(M(indeces_i, indeces_j)).* nu.^(-abs(K)));
    end
end
end


function n = norm_C_to_C(R)
n = sum(abs(R),2);
end

function op_norm = operator_tail_norm(x, nodes, degree)
% TODO : efficiency!
global nu
warning('THIS IS TOOO SLOOOOOW! ! ! ! ! ! ! ! ! ! ')
K_all = -nodes*degree : nodes*degree;
K_center = K_all .* (abs(K_all)<=nodes);
K_tail = K_all .* (abs(K_all)>nodes);

% max |m|>nodes sum |n|<nodes |x_{n-m}| * nu ^(|n|-|m|)
% if n-m>deg*nodes, x_{n-m} = 0

op_norm = 0;
for m = K_tail
    sum_loc = 0;
    for n = K_center
        if abs( n - m ) <=  nodes*degree
            sum_loc = sum_loc + abs( x(n-m+nodes*degree+1)) * nu^(abs(n) - abs(m));
        end
    end
    if sum_loc > op_norm
        op_norm = sum_loc;
    end
end


end