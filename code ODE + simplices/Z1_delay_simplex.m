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
% norm_M = intval(norm_ell1_to_ell1(M_int, x0.size_vector));
% norm_Q = intval(norm_ell1_to_C(Q_int, x0.size_scalar, x0.size_vector));
% first term: A_s * Pi_M ( D_z Pi_infty F)
n_nodes = x0.nodes;
N = alpha0.vector_field.n_equations;

% norm A star 1 with direct computation

nodes_each_section = (degree+1)*x0.nodes;
center_indices = nodes_each_section + 1 + (-x0.nodes:x0.nodes);
tail_indices = nodes_each_section + 1 + [-nodes_each_section: -x0.nodes-1, x0.nodes+1:nodes_each_section];
size_each_section = 2 * nodes_each_section + 1;
K = -nodes_each_section:nodes_each_section;
% SLOW OR MEMORY ISSUES? 
STAR1 = intval(sparse( size_each_section*x0.size_vector, ...
    size_each_section*x0.size_vector));
STAR3 = intval(sparse( (x0.nodes*2+1)*x0.size_vector, ...
    size_each_section*x0.size_vector));
norm_STAR1 = intval(zeros(alpha0.vector_field.n_equations));

nodes = @(vec) (length(vec) -1)/2;
make_conv_mat = @(vec) toeplitz(...
                        [vec(nodes(vec)+1:end);zeros(nodes(vec),1)],...
                        [vec(nodes(vec)+1:-1:1);zeros(nodes(vec),1)]);

for i =1:N
    index_i = (i-1)*size_each_section + (1:size_each_section);
    % index_i_small = (i-1)*(x0.nodes*2+1) + (1:(x0.nodes*2+1));
    % center_indices_i = (i-1)*size_each_section + center_indices;
    center_indices_i_small = (i-1)*(x0.nodes*2+1) + (1:2*x0.nodes+1);
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
                delay_mat = diag(exp(-1i*delay*phi_int*K));
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
norm_A_STAR1 = [norm_ell1_to_ell1(P_int * kron(D3F2_int,ones(1,1+2*nodes_each_section)) * STAR1, x0.size_vector)
    norm_ell1_to_C(R_int * kron(D3F2_int,ones(1,1+2*nodes_each_section)) * STAR1,x0.size_scalar, x0.size_vector)
    norm_STAR1];

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
% norm_STAR3 = intval(zeros(alpha0.vector_field.n_equations));
% for i =1:N
%     for j = 1:N
%         norm_STAR3(i,j) = intval(operator_tail_norm(squeeze(DF.derivative_Fx_toeplix(j,i,:)), ...
%             n_nodes, degree));
%         delay_term_ij = 0;
%         for term = 1:size(DF.derivative_Fx_delay{i,j})
%             if ~isempty(DF.derivative_Fx_delay{i,j}{term})
%                 conv = DF.derivative_Fx_delay{i,j}{term}.convolution;
%                 norm_conv = operator_tail_norm(conv, n_nodes, ...
%                     degree);
%                 delay_term_ij = delay_term_ij +  norm_conv;
%             end
%         end
%         norm_STAR3(i,j) = norm_STAR3(i,j) + delay_term_ij;
%     end
% end
% % I think??
% norm_Kinv_STAR3 = 1/intval(n_nodes+1) * norm_STAR3;

norm_A_STAR3 = [norm_ell1_to_ell1(M_int * STAR3, x0.size_vector)
    norm_ell1_to_C(Q_int * STAR3, x0.size_scalar, x0.size_vector)
    zeros(x0.size_vector,x0.size_vector)];


% STAR4
norm_STAR4 = intval(zeros(alpha0.vector_field.n_equations));
N = alpha0.vector_field.n_equations;
for i =1:N
    for j = 1:N
        norm_STAR4(i,j) = intval(norm(Xi_vector([],DF.derivative_Fx_toeplix(j,i,:))));
        delay_term_ij = 0;
        for term = 1:size(DF.derivative_Fx_delay{i,j})
            if ~isempty(DF.derivative_Fx_delay{i,j}{term})
                conv = DF.derivative_Fx_delay{i,j}{term}.convolution;
                norm_conv = intval(norm(Xi_vector([], conv.')));
                delay_term_ij = delay_term_ij +  norm_conv;
            end
        end
        norm_STAR4(i,j) = norm_STAR4(i,j) + delay_term_ij;
    end
end
norm_Kinv_STAR4 = 1/intval(n_nodes+1) * norm_STAR4;

norm_A_STAR4 = [norm_P * norm_D3F2 * abs(phi_int)^-1 * norm_Kinv_STAR4
    norm_R * norm_D3F2 * abs(phi_int)^-1 * norm_Kinv_STAR4
    abs(phi_int)^-1 * norm_Kinv_STAR4];


Z1 = sum([norm_A_STAR1, norm_A_STAR2, norm_A_STAR3 + norm_A_STAR4],2);
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
% TODO : efficiency!
global nu
warning('THIS IS TOOO SLOOOOOW! ! ! ! ! ! ! ! ! ! ')
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
%     for n = K_center
%         if abs( n - m ) <=  nodes*degree
%             sum_loc = sum_loc + abs( x(n-m+nodes*degree+1)) * nu^(abs(n) - abs(m));
%         end
%     end
    if sum_loc > op_norm
        op_norm = sum_loc;
    end
end


end