function Afullop = fullopA_Hopf(A_center, X)

if isa(A_center,'double')
    A_center = compose_overoperator(A_center, X.size_real, X.node_space*2, X.node_time*2);
elseif ~isa(A_center,'overoperator')
    error('Unexpected input type')
end

l1 = X.parameters(2);
l2 = X.parameters(1);

c4 = -l1*l2;
c2 = l1;

denom_tail11 = polynomial([-l2, 0, 1, 0, 0]);
denom_tail22 = polynomial([c4, 0, c2, 0, 0]);
zero_pol = polynomial(0);
one_pol = polynomial(1);

pol_tail11 = polynomial_fraction(one_pol, denom_tail11);
pol_tail22 = polynomial_fraction(one_pol, denom_tail22);

cell_mat = cell(2);
cell_mat{1,1} = pol_tail11;
cell_mat{1,2} = zero_pol;
cell_mat{2,1} = zero_pol;
cell_mat{2,2} = pol_tail22;

A_tail_mat = tail_matrix(cell_mat);

Afullop = full_operator(A_center, A_tail_mat);

end