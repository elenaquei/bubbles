function Y = Y_delay_refinement(alpha0, alpha1, alpha2, x0, x1, x2, A0_struct, A1_struct, A2_struct, upper_bound)
Y = compute_Y(alpha0, alpha1, alpha2, x0, x1, x2, A0_struct, A1_struct, A2_struct);
depth = 1;
if any(Y > upper_bound)
    Y = refine(alpha0, alpha1, alpha2, x0, x1, x2, A0_struct, A1_struct, A2_struct, upper_bound, depth);
end

end


function Y = refine(alpha0, alpha1, alpha2, x0, x1, x2, A0_struct, A1_struct, A2_struct, upper_bound, depth)
if depth > 1
    fprintf('We need %i refinements\n',depth)
    if depth > 6
        error('Too many refinements')
    end
end
x01 = 0.5 * (x0 + x1);
x02 = 0.5 * (x0 + x2);
x12 = 0.5 * (x2 + x1);

alpha01 = alpha0;
alpha02 = alpha0;
alpha12 = alpha0;

for i = 1:3
    alpha01.scalar_equations.linear_coef{i} = (alpha1.scalar_equations.linear_coef{i} + alpha0.scalar_equations.linear_coef{i})/2;
    alpha02.scalar_equations.linear_coef{i} = (alpha0.scalar_equations.linear_coef{i} + alpha2.scalar_equations.linear_coef{i})/2;
    alpha12.scalar_equations.linear_coef{i} = (alpha1.scalar_equations.linear_coef{i} + alpha2.scalar_equations.linear_coef{i})/2;
end

A01_struct = mean_struct(A0_struct, A1_struct);
A02_struct = mean_struct(A0_struct, A2_struct);
A12_struct = mean_struct(A2_struct, A1_struct);

Y0 = Y_delay_simplex(alpha0, alpha01, alpha02, x0, x01, x02, A0_struct, A01_struct, A02_struct);
if any(Y0>upper_bound)
    Y0 = refine(alpha0, alpha01, alpha02, x0, x01, x02, A0_struct, A01_struct, A02_struct, upper_bound, depth+1);
end
Y1 = Y_delay_simplex(alpha01, alpha1, alpha12, x01, x1, x12, A01_struct, A1_struct, A12_struct);
if any(Y1>upper_bound)
    Y1 = refine(alpha01, alpha1, alpha12, x01, x1, x12, A01_struct, A1_struct, A12_struct, upper_bound, depth+1);
end
Y2 = Y_delay_simplex(alpha02, alpha12, alpha2, x02, x12, x2, A02_struct, A12_struct, A2_struct);
if any(Y2>upper_bound)
    Y2 = refine(alpha02, alpha12, alpha2, x02, x12, x2, A02_struct, A12_struct, A2_struct, upper_bound, depth+1);
end
Y3 = Y_delay_simplex(alpha02, alpha12, alpha01, x02, x12, x01, A02_struct, A12_struct, A01_struct);
if any(Y3>upper_bound)
    Y3 = refine(alpha02, alpha12, alpha01, x02, x12, x01, A02_struct, A12_struct, A01_struct, upper_bound, depth+1);
end

Y = max(max(Y3,Y0), max(Y1, Y2));

end


function A12_struct = mean_struct(A0_struct, A1_struct)
[A0, M0, P0, Q0, R0, phi0, D3F20] = extract_from_struct(A0_struct);
[A1, M1, P1, Q1, R1, phi1, D3F21] = extract_from_struct(A1_struct);
A12_struct.A = (A0+A1)/2;
A12_struct.M = (M0+M1)/2;
A12_struct.P = (P0+P1)/2;
A12_struct.Q = (Q0+Q1)/2;
A12_struct.R = (R0+R1)/2;
A12_struct.phi = (phi0+phi1)/2;
A12_struct.D3F2 = (D3F20+D3F21)/2;
end


function [A, M, P, Q, R, phi, D3F2] = extract_from_struct(A_struct)
A = A_struct.A;
M = A_struct.M;
P = A_struct.P;
Q = A_struct.Q;
R = A_struct.R;
phi = A_struct.phi;
D3F2 = A_struct.D3F2;
end

function Y = compute_Y(alpha0, alpha1, alpha2, x0, x1, x2, A0_struct, A1_struct, A2_struct)
global use_intlab

if use_intlab
    [alpha0, alpha1, alpha2, x0, x1, x2] = everything_intval(alpha0, alpha1, alpha2, x0, x1, x2);
end

alpha  = interpolation(alpha0, alpha1, alpha2);
n_nodes_long = alpha.vector_field.deg_vector*x0.nodes;

x_int = interpolation(x0, x1, x2);

bool_long = 1;

[A0, M0, P0, Q0, R0, phi0, D3F20] = extract_from_struct(A0_struct);
[A1, M1, P1, Q1, R1, phi1, D3F21] = extract_from_struct(A1_struct);
[A2, M2, P2, Q2, R2, phi2, D3F22] = extract_from_struct(A2_struct);
big_A0 = create_A_of_size(A0, M0, P0, Q0, R0, phi0, D3F20, n_nodes_long);
big_A1 = create_A_of_size(A1, M1, P1, Q1, R1, phi1, D3F21, n_nodes_long);
big_A2 = create_A_of_size(A2, M2, P2, Q2, R2, phi2, D3F22, n_nodes_long);
    
big_A_int = interpolation(big_A0, big_A1, big_A2);

    function Y_0 = brute_force_Y(A, alpha, x)
        F = Xi_vec2vec(apply(alpha, x, bool_long));
        Y_0 = norm(vec2Xi_vec(A*F, x.size_scalar, x.size_vector));
    end


Y = brute_force_Y(big_A_int, alpha, x_int);
end