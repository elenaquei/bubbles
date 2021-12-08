function block = block_norm(M, P, Q, R, D3F2, psi, x)
global norm_weight

if isempty(norm_weight)
    norm_weight = ones(x.size_scalar+x.size_vector,1);
elseif length(norm_weight) ~= x.size_scalar+x.size_vector
    if length(norm_weight) == 1
        norm_weight = norm_weight * ones(x.size_scalar+x.size_vector,1);
    else
        error('The weight of the norm is incompatible with the size of the problem')
    end
end

size_scalar = x.size_scalar;
size_vector = x.size_vector;
nodes = x.nodes;

scalar_weights = norm_weight(1:size_scalar);
vector_weights = norm_weight(size_scalar+1:end);

abs_psi = norm_weight(1) * abs(psi) / norm_weight(1);
abs_D3F2 = diag(scalar_weights) * abs(D3F2) * diag(vector_weights.^-1);
norm_P = diag(vector_weights) * norm_C_to_ell1(P, nodes, size_scalar, size_vector) * diag(scalar_weights.^-1);
norm_Q = diag(scalar_weights) * norm_ell1_to_C(Q, nodes, size_scalar, size_vector)* diag(vector_weights.^-1);
norm_M = diag(vector_weights) * norm_ell1_to_ell1(M, nodes, size_vector) * diag(vector_weights.^-1);
norm_R = diag(scalar_weights) * abs(R) * diag(scalar_weights.^-1);
norm_V1 = norm_P * abs_D3F2 / abs_psi;
norm_V2 = norm_R * abs_D3F2 / abs_psi;
norm_V3 = 1/ abs_psi;

block = [norm_R, max(norm_Q,norm_V2);
    norm_P, max(norm_M, norm_V1 + norm_V3)];
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
