function [A_small, M, P, Q, R, phi, D3F2, DF_small] = A_delay_symplex(alpha, x)
% function [A_small, M, P, Q, R, phi, D3F2, DF_small] = A_delay_symplex(alpha, x)
%
% returns the numerical center of the A operator, and all the bits and
% pieces to construct the A operator up to any size
global use_intlab
temp_use_intlab = use_intlab;
use_intlab = 0;

small_der = derivative(alpha, x, 0);
DF_small = derivative_to_matrix(small_der);

inverse_mat_der = inv(DF_small);

index_scalar = 1:x.size_scalar;
index_vector = x.size_scalar+1:length(x);

M = inverse_mat_der(index_vector,index_vector);
P = inverse_mat_der(index_vector,index_scalar);
Q = inverse_mat_der(index_scalar, index_vector);
R = inverse_mat_der(index_scalar, index_scalar);

phi = x.scalar(1);

% attention, attention, attention: D3F2, as defined, is a hell of a mess!
D3F2 = zeros(alpha.scalar_equations.num_equations, x.size_vector);
if isintval(small_der.derivative_G_hat)
    D3F2 = intval(D3F2);
end
index_non_comp_DF = alpha.scalar_equations.number_equations_lin+alpha.scalar_equations.number_equations_pol+1 : alpha.scalar_equations.num_equations;
index_non_comp_DG = alpha.scalar_equations.number_equations_pol+1:alpha.scalar_equations.number_equations_pol+alpha.scalar_equations.number_equations_non_computable;
D3F2(index_non_comp_DF) = small_der.derivative_G_hat(index_non_comp_DG);

A_small = inverse_mat_der;

use_intlab = temp_use_intlab;