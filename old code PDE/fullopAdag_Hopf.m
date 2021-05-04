function Adagfullop = fullopAdag_Hopf(Adag_center, X)

if isa(Adag_center,'double')
    Adag_center = compose_overoperator(Adag_center, X.size_real, X.node_space*2, X.node_time*2);
elseif ~isa(Adag_center,'overoperator')
    error('Unexpected input type')
end

l1 = X.parameters(2);
l2 = X.parameters(1);

c4 = -l1*l2;
c2 = l1;

pol_tail1 = polynomial([-l2, 0, 1, 0, 0]);
pol_tail2 = polynomial([c4, 0, c2, 0, 0]);
zero_pol = polynomial(0);

cell_mat = cell(2);
cell_mat{1,1} = pol_tail1;
cell_mat{1,2} = zero_pol;
cell_mat{2,1} = zero_pol;
cell_mat{2,2} = pol_tail2;

Adag_tail_mat = tail_matrix(cell_mat);

Adagfullop = full_operator(Adag_center, Adag_tail_mat);

end