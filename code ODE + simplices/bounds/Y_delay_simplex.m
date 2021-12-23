function Y = Y_delay_simplex(alpha0, alpha1, alpha2, x0, x1, x2, A0_struct, A1_struct, A2_struct)
global use_intlab

if use_intlab
    [alpha0, alpha1, alpha2, x0, x1, x2] = everything_intval(alpha0, alpha1, alpha2, x0, x1, x2);
end

alpha  = interpolation(alpha0, alpha1, alpha2);
n_nodes_long = alpha.vector_field.deg_vector*x0.nodes;

x_int = interpolation(x0, x1, x2);

bool_long = 1;
Fx_0 = apply(alpha, x0, bool_long);

    function [A, M, P, Q, R, phi, D3F2] = extract_from_struct(A_struct)
        A = A_struct.A;
        M = A_struct.M;
        P = A_struct.P;
        Q = A_struct.Q;
        R = A_struct.R;
        phi = A_struct.phi;
        D3F2 = A_struct.D3F2;
    end

[A0, M0, P0, Q0, R0, phi0, D3F20] = extract_from_struct(A0_struct);
[A1, M1, P1, Q1, R1, phi1, D3F21] = extract_from_struct(A1_struct);
[A2, M2, P2, Q2, R2, phi2, D3F22] = extract_from_struct(A2_struct);
big_A0 = create_A_of_size(A0, M0, P0, Q0, R0, phi0, D3F20, n_nodes_long);
big_A1 = create_A_of_size(A1, M1, P1, Q1, R1, phi1, D3F21, n_nodes_long);
big_A2 = create_A_of_size(A2, M2, P2, Q2, R2, phi2, D3F22, n_nodes_long);
    
    
    function Y_0 = non_cont_Y(A, alpha, x)
        F = Xi_vec2vec(apply(alpha, x, bool_long));
        Y_0 = norm(vec2Xi_vec(A*F, x.size_scalar, x.size_vector));
    end


Y0 = non_cont_Y(big_A0, alpha0, x0);
Y1 = non_cont_Y(big_A1, alpha1, x1);
Y2 = non_cont_Y(big_A2, alpha2, x2);

Y_edge = intval(max(Y0, max(Y1,Y2)));

% first term Taylor series
% first der w.r.t. s
first_s_der = intval(0*reshape(x0,n_nodes_long));
x_int_short = reshape(x_int, x0.nodes);
x_Delta1 = x1 - x0;
x_Delta2 = x2 - x0;
coefDelta1 = cell(3,1);
coefDelta2 = cell(3,1);
coef0 = alpha0.scalar_equations.linear_coef;
coef1 = alpha1.scalar_equations.linear_coef;
coef2 = alpha2.scalar_equations.linear_coef;
coef_int = alpha.scalar_equations.linear_coef;
for i = 1:3
    coefDelta1{i} = coef1{i} - coef0{i};
    coefDelta2{i} = coef2{i} - coef0{i};
    coef_int{i} = interpolation(coef0{i}, coef1{i}, coef2{i});
end
coef_int{3} = 0* coef_int{3};
first_s_der.scalar(1:alpha1.scalar_equations.number_equations_lin) = ...
    interpolation(apply_lin_coef(coefDelta1, x_int_short), apply_lin_coef(coefDelta2, x_int_short)) + ...
    interpolation(apply_lin_coef(coef_int, x_Delta1), apply_lin_coef(coef_int, x_Delta2));

% first der w.r.t. u
DF = derivative(alpha0, x0, bool_long);
DF_mat = derivative_to_matrix(DF, n_nodes_long);
x_Delta1_long = reshape(x_Delta1, n_nodes_long);
x_Delta2_long = reshape(x_Delta2, n_nodes_long);
first_u_der = vec2Xi_vec(interpolation(DF_mat*Xi_vec2vec(x_Delta1_long),...
    DF_mat*Xi_vec2vec(x_Delta2_long)), x0.size_scalar, x0.size_vector, n_nodes_long);

big_A_int = interpolation(big_A0, big_A1, big_A2);

Y_der = big_A_int * Xi_vec2vec(first_s_der + first_u_der);

Y_der_norm = intval(norm(vec2Xi_vec(Y_der, x0.size_scalar, x0.size_vector, n_nodes_long)));


% residual
mixed_second_der = intval(zeros(length(x0),1));
mixed_second_der(1:alpha1.scalar_equations.number_equations_lin) = ...
    interpolation(...
    interpolation(apply_lin_coef(coefDelta1, x_Delta1), apply_lin_coef(coefDelta2, x_Delta2)),...
    interpolation(apply_lin_coef(coefDelta2, x_Delta1), apply_lin_coef(coefDelta1, x_Delta2)));
mixed_second_der = sup(abs(mixed_second_der));

norm_x_Delta = max( sup(norm(x_Delta1)), sup(norm(x_Delta2)));

bound_second_der = DDF(alpha, x_int); 
% TODO: should this be weighted? the norm of XDelta is already weighted...

A_int = interpolation(A0, A1, A2);
M_int = interpolation(M0, M1, M2);
P_int = interpolation(P0, P1, P2);
Q_int = interpolation(Q0, Q1, Q2);
R_int = interpolation(R0, R1, R2);
phi_int = interpolation(phi0, phi1, phi2);
D3F2_int = interpolation(D3F20, D3F21, D3F22);
normA = block_norm(M_int, P_int, Q_int, R_int, D3F2_int/(x0.nodes+1), phi_int, x0); % weighted

norm_AdxdsF = norm(vec2Xi_vec(A_int * mixed_second_der,x0));

residual = norm_AdxdsF + normA *  bound_second_der*max(norm_x_Delta.^2);

Y = sup(Y_edge + Y_der_norm + residual);

end


function y_scal = apply_lin_coef(linear_coef, xi_vec)
global use_intlab
number_equations_lin = length(linear_coef{3});
y_scal = zeros(number_equations_lin,1);
if use_intlab || isintval(xi_vec) || isintval(linear_coef{1})
    y_scal = intval(y_scal);
end

% linear equations
if number_equations_lin>0
    if (isempty(use_intlab) || ~use_intlab) && ~isintval(xi_vec) && ~isintval(linear_coef{1})
        temp=zeros(size(linear_coef{3}));
        for i=1:number_equations_lin
            temp(i)=sum(sum(squeeze(linear_coef{2}(i,:,:)).*xi_vec.vector,1));
        end
        y_scal(1:number_equations_lin)=linear_coef{3}+(linear_coef{1}*(xi_vec.scalar.'))+temp;
    else
        temp=intval(zeros(size(linear_coef{3})));
        for i=1:number_equations_lin
            temp(i)=sum(sum(squeeze(intval(linear_coef{2}(i,:,:))).*xi_vec.vector,1));
        end
        y_scal(1:number_equations_lin)=intval(linear_coef{3})+(intval(linear_coef{1})*(xi_vec.scalar.'))+temp;
    end
end
end



function DDF = DDF(alpha, x)
global nu
global norm_weight

modes = x.nodes;
if nu^(2*modes+2) < exp(1)
    Lemma_bound = nu^(2*modes+2)/ (exp(1) * log(nu^(2*modes+2)));
else
    Lemma_bound = 1;
end

M = x.size_scalar;
N = x.size_vector;
I = eye(M+N);

K = -x.nodes:x.nodes;
x_norm = norm(x).';
lambda_norm = x_norm(1:x.size_scalar);
v_norm = x_norm(x.size_scalar+1:end);
Kx = x;
for j = 1:x.size_vector
    Kx.vector(j,:) = K.*x.vector(j,:);
end
Kx_norm = norm(Kx);
Kv_norm = Kx_norm(x.size_scalar+1:end);
KKx = x;
for j = 1:x.size_vector
    KKx.vector(j,:) = K.*Kx.vector(j,:);
end
KKx_norm = norm(KKx);
KKv_norm = KKx_norm(x.size_scalar+1:end);

DDF = intval(zeros(length(x_norm),1));

for i=1:alpha.vector_field.n_equations % equation
    DpsiDF_i = 0;
    for j=1:alpha.vector_field.n_terms(i) % element of the equation
        power_loc=[alpha.vector_field.power_scalar{i}(:,j).', alpha.vector_field.power_vector{i}{j}.'];
        const=abs(alpha.vector_field.value{i}(j));
        delay_loc = alpha.vector_field.delay{i}{j};
        
        X = prod(x_norm.^power_loc);
        V_delay = prod(v_norm(find(delay_loc)));
        
        if any(delay_loc)
            delay_term1 = 0;
            delay_term2 = 0;
            if power_loc(1)>0
                for l = 1:length(delay_loc)
                    if delay_loc(l) == 0
                        continue
                    end
                    delay_term1 = delay_term1 + abs(delay_loc(l))/v_norm(l) * ( Kv_norm(l) + KKv_norm(l));
                end
                delay_term1 = power_loc(1)/x_norm(1) *delay_term1;
            end
            for l1 = 1:length(delay_loc)
                if delay_loc(l1) == 0
                    continue
                end
                for l2 = 1:length(delay_loc)
                    if delay_loc(l2) == 0
                        continue
                    end
                    if l1 == l2
                        continue
                    end
                    delay_term2 = delay_term2 + abs(delay_loc(l1))/v_norm(l1) *  Kv_norm(l1) * abs(delay_loc(l2))/v_norm(l2) *  Kv_norm(l2) ;
                end
            end
        end
        
        non_delay_term = 0;
        delay_term3 = 0;
        for k = 1:length(x_norm)
            k_is_1 = 1*(k==1);
            non_delay_term = non_delay_term + power_loc(k) /x_norm(k) *(power_loc(1) - k_is_1)/x_norm(1);
            for l = 1:length(delay_loc)
                if delay_loc(l) == 0
                    continue
                end
                delay_term3 = delay_term3 +  power_loc(k) /x_norm(k) * Kv_norm(l) * abs(delay_loc(l))/v_norm(l);
            end
        end
        if all(delay_loc == 0)
            DpsiDF_ij = const * X * V_delay * non_delay_term;
            DpsiDF_i = DpsiDF_i + DpsiDF_ij;
            continue
        end
        delay_term4 = 0;
        for k = 1:length(delay_loc)
            if delay_loc(k) == 0
                continue
            end
            delay_term5 = sum(Kv_norm .* abs(delay_loc)./v_norm.');
            
            delay_term4 = delay_term4 + 1/ v_norm(k) * (power_loc(1)/x_norm(1) + Lemma_bound * abs(delay_loc(k)) + delay_term5);
        end
        DpsiDF_ij = const * X * V_delay * ( delay_term1 + delay_term2 + delay_term3 + delay_term4 + non_delay_term);
        DpsiDF_i = DpsiDF_i + DpsiDF_ij;
        if any(DpsiDF_ij<0)
            error('An absolute value is missing!')
        end
    end
    DDF(i+x.size_scalar) = DpsiDF_i;
end

% second derivative of G(x) = 0

% second derivative of the scalar polynomial equations
for i=1:alpha.scalar_equations.number_equations_pol % equation
    DpsiDF_i = 0;
    for j=1:alpha.scalar_equations.polynomial_equations.n_terms(i) % element of the equation
        d=[alpha.scalar_equations.polynomial_equations.power_scalar{i}(:,j).', alpha.scalar_equations.polynomial_equations.power_vector{i}{j}.'];
        const=alpha.scalar_equations.polynomial_equations.value{i}(j);
        
        X = prod(x_norm.^power_loc.');
        
        non_delay_term = 0;
        for k = 1:length(x_norm)
            non_delay_term = non_delay_term + d(k) /x_norm(k) *(d(1) - k_is_1)/x_norm(1);
        end
        DpsiDF_ij = const * X * V_delay * non_delay_term;
        DpsiDF_i = DpsiDF_i + DpsiDF_ij;
    end
    DDF(i+alpha.scalar_equations.number_equations_lin) = DpsiDF_i;
end

% % % % % % % pasting second section

for i=1:alpha.vector_field.n_equations % equation
    DDF_i = 0;
    for j=1:alpha.vector_field.n_terms(i) % element of the equation
        power_loc=[alpha.vector_field.power_scalar{i}(:,j).', alpha.vector_field.power_vector{i}{j}.'];
        const=abs(alpha.vector_field.value{i}(j));
        abs_delay_loc = abs(alpha.vector_field.delay{i}{j}).';
        n_delays = length(find(abs_delay_loc));
        X = prod(x_norm.^power_loc);
        V_delay = prod(v_norm(find(abs_delay_loc)));
        
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
        
        DDF_ij = const * X * V_delay * ( delay_term1 + delay_term2 + delay_term3 + no_delay_term);
        DDF_i = DDF_i + DDF_ij;
    end
    DDF(i+x.size_scalar) = DDF(i+x.size_scalar) + DDF_i;
end

% second derivative of G(x) = 0

% second derivative of the scalar polynomial equations
for i=1:alpha.scalar_equations.number_equations_pol % equation
    % alpha.scalar_equations.number_equations_lin+
    DDF_i = 0;
    for j=1:alpha.scalar_equations.polynomial_equations.n_terms(i) % element of the equation
        power_loc=[alpha.scalar_equations.polynomial_equations.power_scalar{i}(:,j).', alpha.scalar_equations.polynomial_equations.power_vector{i}{j}.'];
        const=alpha.scalar_equations.polynomial_equations.value{i}(j);
        
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
    DDF(i+alpha.scalar_equations.number_equations_lin) = DDF(i+alpha.scalar_equations.number_equations_lin) + DDF_i;
end

% added user input derivative
[~,~,~,DD_user] = alpha.scalar_equations.non_computable_function(x_norm);
index_non_comp = alpha.scalar_equations.num_equations - alpha.scalar_equations.number_equations_non_computable +1 : alpha.scalar_equations.num_equations;
DDF(index_non_comp) = ...
    DDF(index_non_comp) + sum(DD_user,2);

DDF = DDF .* norm_weight;
end


