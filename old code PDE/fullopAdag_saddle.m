function A_fullop = fullopAdag_saddle(A,X)

% A_op = cell(3);
pol_tail = cell(6);

x = cell(3,1);
% [x{1},x{2},x{3}] = split(X);

l1 = X.parameters(1);
l2 = X.parameters(2);

c4 = -l1*l2;
c2 = l1;

zero_pol = polynomial(0);

pol_tail_diag1 = polynomial([-l2, 0, 1, 0, 0]);
pol_tail_diag2 = polynomial([c4, 0, c2, 0, 0]);
pol_tail_else = zero_pol;

% index_1 =1:length(x{1});
for i = 1:3
    for j = 1:3
%         A_tempij = A(length(x{1})*(i-1)+index_1,length(x{1})*(j-1)+index_1);
        if i==j
%             A_split_ij =  big_Adagger(A_tempij,x{i});
            pol_tail{2*(i-1)+1,2*(j-1)+1} = pol_tail_diag1;
            pol_tail{2*(i-1)+2,2*(j-1)+2} = pol_tail_diag2;
            
            pol_tail{2*(i-1)+1,2*(j-1)+2} = pol_tail_else;
            pol_tail{2*(i-1)+2,2*(j-1)+1} = pol_tail_else;
        else
%             A_split_ij = [0*A_tempij, 0*A_tempij, 0*A_tempij
%                 0*A_tempij, A_tempij, 0*A_tempij
%                 0*A_tempij, 0*A_tempij, 0*A_tempij];
            
            pol_tail{2*(i-1)+1,2*(j-1)+1} = pol_tail_else;
            pol_tail{2*(i-1)+2,2*(j-1)+2} = pol_tail_else;
            
            pol_tail{2*(i-1)+1,2*(j-1)+2} = pol_tail_else;
            pol_tail{2*(i-1)+2,2*(j-1)+1} = pol_tail_else;
        end
%         A_op{i,j} = compose_overoperator( A_split_ij, x{1});
    end
end

tail_A = tail_matrix(pol_tail);

center_A= compose_overoperator3(A, X.size_real, X.node_space*2, X.node_time*2);

% center_A = overoperator3(A_op{1,1},A_op{1,2},A_op{1,3},...
%     A_op{2,1},A_op{2,2},A_op{2,3},...
%     A_op{3,1},A_op{3,2},A_op{3,3});

A_fullop = full_operator(center_A,tail_A);
