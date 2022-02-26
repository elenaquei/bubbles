function A_fullop = fullopA_saddle(A,X)


pol_tail = cell(6);

l1 = X.parameters(1);
l2 = X.parameters(2);

c4 = -l1*l2;
c2 = l1;

denom_tail1 = polynomial([-l2, 0, 1, 0, 0]);
denom_tail2 = polynomial([c4, 0, c2, 0, 0]);
zero_pol = polynomial(0);
one_pol = polynomial(1);

pol_tail_diag1 = polynomial_fraction(one_pol, denom_tail1);
pol_tail_diag2 = polynomial_fraction(one_pol, denom_tail2);
pol_tail_else = zero_pol;

for i = 1:3
    for j = 1:3
        if i==j
            pol_tail{2*(i-1)+1,2*(j-1)+1} = pol_tail_diag1;
            pol_tail{2*(i-1)+2,2*(j-1)+2} = pol_tail_diag2;
            
            pol_tail{2*(i-1)+1,2*(j-1)+2} = pol_tail_else;
            pol_tail{2*(i-1)+2,2*(j-1)+1} = pol_tail_else;
        else
            
            pol_tail{2*(i-1)+1,2*(j-1)+1} = pol_tail_else;
            pol_tail{2*(i-1)+2,2*(j-1)+2} = pol_tail_else;
            
            pol_tail{2*(i-1)+1,2*(j-1)+2} = pol_tail_else;
            pol_tail{2*(i-1)+2,2*(j-1)+1} = pol_tail_else;
        end
    end
end

tail_A = tail_matrix(pol_tail);

center_A= compose_overoperator3(A, X.size_real, X.node_space*2, X.node_time*2);

A_fullop = full_operator(center_A,tail_A);
