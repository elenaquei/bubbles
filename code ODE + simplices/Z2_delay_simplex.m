function [Z2]=Z2_delay_simplex(x0,x1,x2,alpha, Rmax)

modes = x0.n_nodes;

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

block_norm_A = block_norm(M_int, P_int, Q_int, R_int, D3F2_int/(modes+1), phi_int, x0);

upper_bound_psi_DF = dpsiDF(alpha, x, Rmax);

Z22 = block_norm_A * upper_bound_psi_DF;

block_norm_AM = block_norm((modes+1)*M_int, P_int, (modes+1)*Q_int, R_int, D3F2_int, phi_int, x0);
DDF_without_dpsi = upper_bound_DDF_without_psi(alpha, x, Rmax);
Z21 = block_norm_AM * DDF_without_dpsi;

Z2 = Z22 + Z21;
end

function DpsiDF = dpsiDF(alpha, x, Rmax)
% MOSTLY DONE
modes = x.nodes;
if nu^(2*modes+2) < exp(1)
    Lemma_bound = nu^(2*modes+2)/ (exp(1) * log(nu^(2*modes+2)));
else
    Lemma_bound = 1;
end

K = -x.nodes:x.nodes;
x_norm = norm(x) + Rmax;
lambda_norm = x_norm(1:x.size_scalar);
v_norm = x_norm(x.size_scalar:end);
Kx = x;
for j = 1:x.size_vector
    Kx.vector(j,:) = K.*x.vector;
end
Kx_norm = norm(Kx)+Rmax;
Kv_norm = Kx_norm(x.size_scalar+1:end);
KKx = x;
for j = 1:x.size_vector
    KKx.vector(j,:) = K.*Kx.vector;
end
KKx_norm = norm(KKx)+Rmax;
KKv_norm = KKx_norm(x.size_scalar+1:end);

DpsiDF = zeros(length(x),1);

for i=1:alpha.vector_field.n_equations % equation
    DpsiDF_i = 0;
    for j=1:alpha.vector_field.n_terms(i) % element of the equation
        power_loc=[alpha.vector_field.power_scalar{i}(:,j).', alpha.vector_field.power_vector{i}{j}.'];
        const=abs(alpha.vector_field.value{i}(j));
        delay_loc = alpha.vector_field.delay{i}{j};
        
        X = prod(x_norm.^power_loc);
        V_delay = prod(v_norm(delay_loc));
        
        if any(delay_loc)
            if power_loc(1)>0
                delay_term1 = 0;
                for l = find(delay_loc)
                    delay_term1 = delay_term1 + abs(delay_loc(l))/v_norm(l) * ( Kv_norm(l) + KKv_norm(l));
                end
                delay_term1 = power_loc(1)/x_norm(1) *delay_term1;
                delay_term2 = 0;
            end
            for l1 = find(delay_loc)
                for l2 = find(delay_loc)
                    if l1 == l2
                        continue
                    end
                    delay_term2 = delay_term2 + abs(delay_loc(l1))/v_norm(l1) *  Kv_norm(l1) * abs(delay_loc(l2))/v_norm(12) *  Kv_norm(l2) ;
                end
            end
        end
        
        non_delay_term = 0;
        delay_term3 = 0;
        for k = 1:length(x)
            non_delay_term = non_delay_term + d(k) /x_norm(k) *(d(1) - k_is_1)/x_norm(1);
            for l = find(delay_loc)
                delay_term3 = delay_term3 +  d(k) /x_norm(k) * Kv_norm(l) * abs(delay_loc(l))/v_norm(l);
            end
        end
        if all(delay_loc == 0)
            DpsiDF_ij = const * X * V_delay * non_delay_term;
            DpsiDF_i = DpsiDF_i + DpsiDF_ij;
            continue
        end
        delay_term4 = 0;
        for k = find(delay_loc)
            delay_term5 = sum(Kv_norm .* delay_loc./v_norm);
            
            delay_term4 = delay_term4 + 1/ v_norm(k) * (power_loc(1)/x_norm(1) + Lemma_bound * delay_loc(k) + delay_term5);
        end
        DpsiDF_ij = const * X * V_delay * ( delay_term1 + delay_term2 + delay_term3 + delay_term4 + non_delay_term);
        DpsiDF_i = DpsiDF_i + DpsiDF_ij;
    end
    DpsiDF(i+x.size_scalar) = DpsiDF_i;
end

% second derivative of G(x) = 0

% second derivative of the scalar polynomial equations
for i=1:alpha.scalar_equations.number_equations_pol % equation
    % alpha.scalar_equations.number_equations_lin+
    DpsiDF_i = 0;
    for j=1:alpha.scalar_equations.polynomial_equations.n_terms(i) % element of the equation
        d=[alpha.scalar_equations.polynomial_equations.power_scalar{i}(:,j).', alpha.scalar_equations.polynomial_equations.power_vector{i}{j}.'];
        const=alpha.scalar_equations.polynomial_equations.value{i}(j);
        N=length(d);
        e=eye(N);
        
        X = prod(x_norm.^power_loc);
        
        non_delay_term = 0;
        for k = 1:length(x)
            non_delay_term = non_delay_term + d(k) /x_norm(k) *(d(1) - k_is_1)/x_norm(1);
        end
        DpsiDF_ij = const * X * V_delay * non_delay_term;
        DpsiDF_i = DpsiDF_i + DpsiDF_ij;
    end
    DpsiDF(i+alpha.scalar_equations.num_equations_lin) = DpsiDF_i;
end

% MISSING: USER INPUT DERIVATIVE!!
warning('User input derivative not included yet')
end

function upper_bound_DDF_without_psi(alpha, x, Rmax)

modes = x.nodes;
if nu^(2*modes+2) < exp(1)
    Lemma_bound = nu^(2*modes+2)/ (exp(1) * log(nu^(2*modes+2)));
else
    Lemma_bound = 1;
end
M = x.size_scalar;
N = x.size_vector;
x_norm = norm(x) + Rmax;
lambda_norm = x_norm(1:x.size_scalar);
v_norm = x_norm(x.size_scalar:end);
I = eye(M+N);
DDF = zeros(length(x),1);

for i=1:alpha.vector_field.n_equations % equation
    DDF_i = 0;
    for j=1:alpha.vector_field.n_terms(i) % element of the equation
        power_loc=[alpha.vector_field.power_scalar{i}(:,j).', alpha.vector_field.power_vector{i}{j}.'];
        const=abs(alpha.vector_field.value{i}(j));
        abs_delay_loc = abs(alpha.vector_field.delay{i}{j});
        n_delays = length(find(abs_delay_loc));
        X = prod(x_norm.^power_loc);
        V_delay = prod(v_norm(abs_delay_loc));
        
        no_delay_term = 0;
        for k = 1:M+N
            no_delay_term = no_delay_term + power_loc(k) / x_norm(k) * sum(1./x_norm.*(power_loc - I(k,:)));
        end
        if all(abs_delay_loc == 0)
            DDF_ij = const * X * no_delay_term;
            DDF_i = DDF_i + DDF_ij;
            continue
        end
        
        delay_term1 = Lemma_bound * sum(abs_delay_loc ./ v_norm);
        delay_term2 = n_delays *(M+N) + ( n_delays + Lemma_bound * sum(abs_delay_loc .* v_norm)) * sum(power_loc);
        delay_term3 = (Lemma_bound * sum(v_norm.*abs_delay_loc) + n_delays) * n_delays;
        
        DDF_ij = const * X * V_delay * ( delay_term1 + delay_term2 + delay_term3 + non_delay_term);
        DDF_i = DDF_i + DDF_ij;
    end
    DDF(i+x.size_scalar) = DDF_i;
end

% second derivative of G(x) = 0

% second derivative of the scalar polynomial equations
for i=1:alpha.scalar_equations.number_equations_pol % equation
    % alpha.scalar_equations.number_equations_lin+
    DpsiDF_i = 0;
    for j=1:alpha.scalar_equations.polynomial_equations.n_terms(i) % element of the equation
        d=[alpha.scalar_equations.polynomial_equations.power_scalar{i}(:,j).', alpha.scalar_equations.polynomial_equations.power_vector{i}{j}.'];
        const=alpha.scalar_equations.polynomial_equations.value{i}(j);
        N=length(d);
        I=eye(N);
        
        X = prod(x_norm.^power_loc);
        
        no_delay_term = 0;
        for k = 1:M+N
            no_delay_term = no_delay_term + power_loc(k) / x_norm(k) * sum(1./x_norm.*(power_loc - I(k,:)));
        end
        if all(abs_delay_loc == 0)
            DDF_ij = const * X * no_delay_term;
            DDF_i = DDF_i + DDF_ij;
            continue
        end
    end
    DDF(i+alpha.scalar_equations.num_equations_lin) = DDF_i;
end

% MISSING: USER INPUT DERIVATIVE!!
warning('User input derivative not included yet')
end

function block = block_norm(M, P, Q, R, D3F2, psi, x)
size_scalar = x.size_scalar;
size_vector = x.size_vector;
nodes = x.nodes;

norm_P = norm_C_to_ell1(P, nodes, size_scalar, size_vector);
norm_Q = norm_ell1_to_C(Q, nodes, size_scalar, size_vector);
norm_M = norm_ell1_to_ell1(M, nodes, size_vector);
norm_R = abs(R);
norm_V1 = norm_P * abs(D3F2) / abs(psi);
norm_V2 = norm_R * abs(D3F2) / abs(psi);
norm_V3 = 1 / abs(psi);

block = [ norm_M + norm_V1 + norm_V3, norm_P
    norm_Q+norm_V2, n_R];
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