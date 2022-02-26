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
inverse_mat_der= symmetrise_A(inverse_mat_der, x);

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
end

function A_sim=symmetrise_A(A,x)
A_sim=A;

length_nodes=2*x.nodes+1;
sym_vec=@(v) ( v + conj(v(end:-1:1,:)))/2;
sym_mat=@(A) ( A + conj(A(end:-1:1,end:-1:1)))/2;

for j=1:x.size_vector
    vec=x.size_scalar+(j-1)*length_nodes+1:x.size_scalar+j*length_nodes;
    all_scal=1:x.size_scalar;
    A_sim(vec,all_scal)=sym_vec(A_sim(vec,all_scal));
end

for jj=1:x.size_vector
    for kk=1:x.size_vector
        vec1=x.size_scalar+(jj-1)*length_nodes+1:x.size_scalar+jj*length_nodes;
        vec2=x.size_scalar+(kk-1)*length_nodes+1:x.size_scalar+kk*length_nodes;
        A_sim(vec1,vec2)=sym_mat(A_sim(vec1,vec2));
    end
end

end