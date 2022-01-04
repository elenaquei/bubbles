function Z1 = Z1_delay_simplex(alpha0, alpha1, alpha2, x0, x1, x2, A0_struct, A1_struct, A2_struct)
global use_intlab
global nu
global norm_weight

if isempty(norm_weight)
    norm_weight = ones(x0.size_scalar+x0.size_vector,1);
elseif length(norm_weight) ~= x0.size_scalar+x0.size_vector
    if length(norm_weight) == 1
        norm_weight = norm_weight * ones(x0.size_scalar+x0.size_vector,1);
    else
        error('The weight of the norm is incompatible with the size of the problem')
    end
end
weight_scalar = norm_weight(1:x0.size_scalar);
weight_vector = norm_weight(x0.size_scalar+1:end);
norm_weight_long = [weight_vector; weight_scalar; weight_vector];

use_intlab = 1;


% we compute the derivative BEFORE we add nodes, so to have an easier
% computation
bool_long = 1;
% we already have extended x to be long enough, we don't need additional
% terms
x_int = interpolation(x0, x1, x2);
alpha_int = interpolation(alpha0, alpha1, alpha2);
DF = derivative(alpha_int, x_int, bool_long);

% reshape all inputs
if x0.nodes<128
    new_nodes = 128;
    % [A_small0, M0, P0, Q0, R0, phi0, D3F20] = create_A_of_size(A0_struct, new_nodes);
    % [A_small1, M1, P1, Q1, R1, phi1, D3F21] = create_A_of_size(A1_struct, new_nodes);
    % [A_small2, M2, P2, Q2, R2, phi2, D3F22] = create_A_of_size(A2_struct, new_nodes);
    
    [~, M0, P0, Q0, R0, phi0, D3F20] = create_A_of_size(A0_struct, new_nodes);
    [~, M1, P1, Q1, R1, phi1, D3F21] = create_A_of_size(A1_struct, new_nodes);
    [~, M2, P2, Q2, R2, phi2, D3F22] = create_A_of_size(A2_struct, new_nodes);
    x0 = reshape(x0, new_nodes);
    % x1 = reshape(x1, new_nodes);
    % x2 = reshape(x2, new_nodes);
    alpha0 = reshape(alpha0, new_nodes);
    % alpha1 = reshape(alpha1, new_nodes);
    % alpha2 = reshape(alpha2, new_nodes);
    DF = reshape(DF, new_nodes);
end

degree = alpha0.vector_field.deg_vector;
n_nodes = x0.nodes;
size_vector = x0.size_vector;
size_scalar = x0.size_scalar;
N = alpha0.vector_field.n_equations;
% n_nodes_long = n_nodes*degree;

% x_int = interpolation(x0, x1, x2);


    function [A, M, P, Q, R, phi, D3F2] = extract_from_struct(A_struct)
        A = A_struct.A;
        M = A_struct.M;
        P = A_struct.P;
        Q = A_struct.Q;
        R = A_struct.R;
        phi = A_struct.phi;
        D3F2 = A_struct.D3F2;
    end

% [A_small0, M0, P0, Q0, R0, phi0, D3F20] = extract_from_struct(A0_struct);
% [A_small1, M1, P1, Q1, R1, phi1, D3F21] = extract_from_struct(A1_struct);
% [A_small2, M2, P2, Q2, R2, phi2, D3F22] = extract_from_struct(A2_struct);
    

% [A_small0, M0, P0, Q0, R0, phi0, D3F20] = A_delay_symplex(alpha0, x0);
% [A_small1, M1, P1, Q1, R1, phi1, D3F21] = A_delay_symplex(alpha1, x1);
% [A_small2, M2, P2, Q2, R2, phi2, D3F22] = A_delay_symplex(alpha2, x2);

% A_small_int = interpolation(A_small0, A_small1, A_small2);

M_int = interpolation(M0, M1, M2);
P_int = interpolation(P0, P1, P2);
Q_int = interpolation(Q0, Q1, Q2);
R_int = interpolation(R0, R1, R2);
phi_int = interpolation(phi0, phi1, phi2);
D3F2_int = interpolation(D3F20, D3F21, D3F22);

norm_D3F2 =diag(weight_scalar) * intval(sup(abs(D3F2_int))) * diag(weight_vector.^-1);
norm_P = diag(weight_vector) * intval(norm_C_to_ell1(P_int, size_scalar, size_vector)) * diag(weight_scalar.^-1);
norm_R = diag(weight_scalar) * abs(R_int) * diag(weight_scalar.^-1);
% norm_M = intval(norm_ell1_to_ell1(M_int, size_vector));
% norm_Q = intval(norm_ell1_to_C(Q_int, size_scalar, size_vector));
% first term: A_s * Pi_M ( D_z Pi_infty F)

% norm A star 1 with direct computation

nodes_each_section = (degree+1)*n_nodes;
center_indices = nodes_each_section + 1 + (-n_nodes:n_nodes);
tail_indices = nodes_each_section + 1 + [-nodes_each_section: -n_nodes-1, n_nodes+1:nodes_each_section];
size_each_section = 2 * nodes_each_section + 1;
K = -nodes_each_section:nodes_each_section;

STAR1 = intval(sparse( size_each_section*size_vector, ...
    size_each_section*size_vector));
STAR3 = intval(sparse( (n_nodes*2+1)*size_vector, ...
    size_each_section*size_vector));
norm_STAR1 = intval(zeros(alpha0.vector_field.n_equations));

nodes = @(vec) (length(vec) -1)/2;
make_conv_mat = @(vec) toeplitz(...
                        [vec(nodes(vec)+1:end);zeros(nodes(vec),1)],...
                        [vec(nodes(vec)+1:-1:1);zeros(nodes(vec),1)]);

for i =1:N
    index_i = (i-1)*size_each_section + (1:size_each_section);
    index_i_small = (i-1)*(n_nodes*2+1) + (1:(n_nodes*2+1));
    % center_indices_i = (i-1)*size_each_section + center_indices;
    center_indices_i_small = (i-1)*(n_nodes*2+1) + (1:2*n_nodes+1);
    tail_indices_i = (i-1)*size_each_section + tail_indices;
    for j = 1:N
        index_j = (j-1)*size_each_section + (1:size_each_section);
        center_indices_j = (j-1)*size_each_section + center_indices;
        tail_indices_j = (j-1)*size_each_section + tail_indices;
        loc_der = squeeze(DF.derivative_Fx_toeplix(j,i,:));
        if length(loc_der)< length(K)
            % padding
            diff = (length(K) - length(loc_der))/2;
            loc_der = [zeros(diff,1); loc_der; zeros(diff,1)];
        elseif length(loc_der) > length(K)
            error('this should never happen')
        end
        convolution_term_ij = make_conv_mat(loc_der);
        delay_term_ij = 0;
        for term = 1:size(DF.derivative_Fx_delay{i,j})
            if ~isempty(DF.derivative_Fx_delay{i,j}{term})
                conv = DF.derivative_Fx_delay{i,j}{term}.convolution;
                if length(conv)< length(K)
                    % padding
                    diff = (length(K) - length(conv))/2;
                    conv = [zeros(diff,1); conv; zeros(diff,1)];
                elseif length(conv) > length(K)
                    error('this should never happen')
                end
                delay = DF.derivative_Fx_delay{i,j}{term}.exp_coef;
                delay_mat = diag(exp(delay*K));
                delay_term_ij = delay_term_ij +  make_conv_mat(conv)*delay_mat;
            end
        end
        convolution_term_ij = (convolution_term_ij + delay_term_ij);
        convolution_term_K_ij = diag(K.^-1) * convolution_term_ij;
        STAR1(tail_indices_i,center_indices_j) = 1i * phi_int^-1 * convolution_term_K_ij(tail_indices,center_indices);
        STAR3(center_indices_i_small,tail_indices_j) = convolution_term_ij(center_indices,tail_indices);
        norm_STAR1(i,j) = norm_ell1_to_ell1(STAR1(index_i, index_j),1);
        % norm_STAR3(i,j) = norm_ell1_to_ell1(STAR3(index_i_small, index_j),1);
    end
end
% clear delay_mat delay_term_ij convolution_term_ij
norm_A_STAR1 = diag(norm_weight_long) * [norm_ell1_to_ell1(P_int * kron(D3F2_int,ones(1,1+2*nodes_each_section)) * STAR1, size_vector)
    norm_ell1_to_C(R_int * kron(D3F2_int,ones(1,1+2*nodes_each_section)) * STAR1,size_scalar, size_vector)
    norm_STAR1] * diag(weight_vector.^-1);

% STAR2
% Pi_infty_indices = ...
 %    [-n_nodes*degree:-n_nodes-1, 0*(-n_nodes:n_nodes), n_nodes+1:n_nodes*degree];
center_indices = n_nodes*degree-1:n_nodes*degree + 2*n_nodes-1;
derivative_parameters_F_tail = DF.derivative_Flambda;
derivative_parameters_F_tail(:,:, center_indices) = 0;

norm_Kinv_STAR2 = intval(zeros(size_vector, size_scalar));
nodes_loc = (size(derivative_parameters_F_tail,3)-1)/2;
Kinv_STAR2 = intval(zeros((2*nodes_loc+1)*size_vector, size_scalar));
norm_Fourier_sequence=@(vec)sup(sum(abs(vec).*(nu.^abs(-(length(vec)-1)/2:(length(vec)-1)/2))));
for i = 1:size_scalar
    for j = 1:size_vector
        indices = (2*nodes_loc+1)*(j-1) + (1:1+2*nodes_loc);
        K_without_zero = -nodes_loc:nodes_loc;
        K_without_zero(K_without_zero==0) = 1;
        K_der_F_tail = (K_without_zero.^-1).*squeeze(derivative_parameters_F_tail(i,j,:)).';
        Kinv_STAR2(indices,i) = K_der_F_tail;
        norm_Kinv_STAR2(j,i) = intval(norm_Fourier_sequence(K_der_F_tail));
    end
end
%norm_Kinv_STAR2 = sum(abs(1./intval(Pi_infty_indices(Pi_infty_indices~=0)))) * norm_STAR2;

norm_A_STAR2 = diag(norm_weight_long) * [norm_C_to_ell1(P_int * kron(D3F2_int,ones(1,1+2*nodes_loc)) * Kinv_STAR2, size_scalar, size_vector)
    abs(R_int * kron(D3F2_int,ones(1,1+2*nodes_loc)) * Kinv_STAR2)
    abs(phi_int)^-1 * norm_Kinv_STAR2] * diag(weight_scalar.^-1);
%[norm_P * norm_D3F2 * abs(phi_int)^-1 * norm_Kinv_STAR2
%    norm_R.' * norm_D3F2 * abs(phi_int)^-1 * norm_Kinv_STAR2
%    abs(phi_int)^-1 * norm_Kinv_STAR2];

% STAR 3
norm_A_STAR3 = diag(norm_weight_long) * [norm_ell1_to_ell1(M_int * STAR3, size_vector)
    norm_ell1_to_C(Q_int * STAR3, size_scalar, size_vector)
    zeros(size_vector,size_vector)] * diag(weight_vector.^-1);


% STAR4
norm_STAR4 = intval(zeros(alpha0.vector_field.n_equations));
N = alpha0.vector_field.n_equations;

for i =1:N
    for j = 1:N
        norm_STAR4_ij = intval(norm_Fourier_sequence(squeeze(DF.derivative_Fx_toeplix(j,i,:)).'));
        delay_term_ij = 0;
        for term = 1:size(DF.derivative_Fx_delay{i,j})
            if ~isempty(DF.derivative_Fx_delay{i,j}{term})
                conv = DF.derivative_Fx_delay{i,j}{term}.convolution;
                norm_conv = intval(norm_Fourier_sequence(conv.'));
                delay_term_ij = delay_term_ij +  norm_conv;
            end
        end
        norm_STAR4(i,j) = norm_STAR4_ij + delay_term_ij;
    end
end
norm_Kinv_STAR4 = 1/intval(n_nodes+1) * ...
    diag(weight_vector) * norm_STAR4 *  diag(weight_vector.^-1);

norm_A_STAR4 = [norm_P * norm_D3F2 * abs(phi_int)^-1 * norm_Kinv_STAR4
    norm_R * norm_D3F2 * abs(phi_int)^-1 * norm_Kinv_STAR4
    abs(phi_int)^-1 * norm_Kinv_STAR4];

Z1_mat = [norm_A_STAR1, norm_A_STAR2, norm_A_STAR3 + norm_A_STAR4];
Z1_finite_vector = Z1_mat(1:size_vector,1:size_vector) + Z1_mat(size_vector+size_scalar+1:end, 1:size_vector);
Z1_infinite_vector = Z1_mat(1:size_vector,size_vector+size_scalar+1:end) + Z1_mat(size_vector+size_scalar+1:end,size_vector+size_scalar+1:end);
Z1_vector = sum(max(Z1_finite_vector,Z1_infinite_vector),2);

Z1_components = sum([norm_A_STAR1, norm_A_STAR2, norm_A_STAR3 + norm_A_STAR4],2);
Z1_scalar = Z1_components(size_vector+(1:size_scalar));
%Z1_vector = Z1_components(1:size_vector) + Z1_components(size_vector+size_scalar+1:end);
Z1 = sup([Z1_scalar;Z1_vector]);
end

function n = norm_C_to_ell1(P, size_scalar, size_vector)
global use_intval
global nu
n = zeros(size_vector,size_scalar);
if isintval(P) || use_intval
    n = intval(n);
end
nodes = (size(P,1)/size_vector - 1)/2;
K = (-nodes:nodes).';
K_big = repmat(K, 1, size_scalar);
for i = 1:size_vector
    indeces = (1:2*nodes+1) + (i-1) * (2*nodes+1);
    n(i,:) = sum(abs(P(indeces, :)).* nu.^(abs(K_big)),1);
end
end

function n = norm_ell1_to_C(Q, size_scalar, size_vector)
global use_intval
global nu
n = zeros(size_scalar, size_vector);
if isintval(Q) || use_intval
    n = intval(n);
end
nodes = (size(Q,2)/size_vector-1)/2;
K = (-nodes:nodes);
K_big = repmat(K, size_scalar, 1);
for i = 1:size_vector
    indeces = (1:2*nodes+1) + (i-1) * (2*nodes+1);
     n(:,i) = max(abs(Q(:, indeces)).* nu.^(-abs(K_big)), [], 2);
end
end

function n = norm_ell1_to_ell1(M, size_vector)
global use_intval
global nu
n = zeros(size_vector, size_vector);
if isintval(M) || use_intval
    n = intval(n);
end
nodes_1 = (size(M,1)/size_vector - 1)/2;
nodes_2 = (size(M,2)/size_vector - 1)/2;
K = (-nodes_1:nodes_1);
J = (-nodes_2:nodes_2);
for i = 1:size_vector
    indeces_i = (1:2*nodes_1+1) + (i-1) * (2*nodes_1+1);
    for j = 1:size_vector
        indeces_j = (1:2*nodes_2+1) + (j-1) * (2*nodes_2+1);
        n(i,j) = max(nu.^(abs(K))*abs(M(indeces_i, indeces_j)).* nu.^(-abs(J)));
    end
end
end


function n = norm_C_to_C(R)
n = sum(abs(R),2);
end

function op_norm = operator_tail_norm(x, nodes, degree)
global nu
K_all = -nodes*degree : nodes*degree;
K_center = K_all(abs(K_all)<=nodes);
K_tail = K_all(abs(K_all)>nodes);

% max |m|>nodes sum |n|<nodes |x_{n-m}| * nu ^(|n|-|m|)
% if n-m>deg*nodes, x_{n-m} = 0

op_norm = 0;
for m = K_tail
    min_index = max(K_all(1),-m-nodes);
    max_index = min(K_all(end),-m+nodes);
    bad_indices_neg = -m-nodes : K_all(1)-1;
    bad_indices_pos = K_all(end)+1:-m+nodes;
    if ~isempty(bad_indices_neg)
        K_center_loc = K_center(length(bad_indices_neg)+1:end);
    elseif ~isempty(bad_indices_pos)
        K_center_loc = K_center(1:end-length(bad_indices_pos));
    else
        K_center_loc = K_center;
    end
    sum_loc = sum( abs(x((min_index:max_index)+nodes*degree+1)).*nu.^(abs(K_center_loc.') - abs(m)));
    if sum_loc > op_norm
        op_norm = sum_loc;
    end
end


end