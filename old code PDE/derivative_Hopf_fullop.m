function DHopf = derivative_Hopf_fullop(X, G, Es)
% function DHopf = derivative_Hopf_fullop(x_PDE, G, Es)
%
% INPUT 
% x_PDE     PDE_vector with 3 unknown parameters, one 1D vector and one
%           2D vector
% Es        PDE_vector or vector compatible with the previous PDE_vector,
%           continuation equation
% G         matrix of 2 lines, each line compatible with the
%           PDE_vector
% Es        (DEFAULT empty) elements for the continuation equation
% flagFullSize (DEFAULT 0) returns as big as possible
% OUTPUT
% DHopf     matrix, derivative of Kuramoto-shivashinki Hopf in x_PDE
% 

DHopf_op = derivative_Hopf(X, G, Es, 1,1);

l1 = X.parameters(1);
l2 = X.parameters(2);
a = X.parameters(3);

y = X.vector1D;
z = squeeze(X.vector2D);

y_in_z = embedding(y,z);
selected_indeces = X.node_time+1+(0:2*X.node_space)*(2*X.node_time+1);
z_in_y = z(selected_indeces);

c4 = -l1*l2;
c2 = l1;

pol_tail_F1 = polynomial([-l2, 0, 1, 0, 0]);
pol_tail_F2 = polynomial([c4, 0, c2, 0, 0]);
zero_pol = polynomial(0);

cell_lin = cell(2);
cell_lin{1,1} = pol_tail_F1;
cell_lin{1,2} = zero_pol;
cell_lin{2,1} = zero_pol;
cell_lin{2,2} = pol_tail_F2;

pol_nonlin = polynomial([1,0]);

cell_nonlin = cell(2);
cell_nonlin{1,1} = pol_nonlin;
cell_nonlin{1,2} = zero_pol;
cell_nonlin{2,1} = pol_nonlin;
cell_nonlin{2,2} = pol_nonlin;

coef_nonlin = cell(2);
coef_nonlin{1,1} = 2*1i*l1*y;
coef_nonlin{1,2} = 0;
coef_nonlin{2,1} = 4*1i*z_in_y;
coef_nonlin{2,2} = 4*1i*(y_in_z + 2*a*z);% check!

DHopf_tail_lin = tail_matrix(cell_lin);
DHopf_tail_nonlin = tail_matrix(cell_nonlin);

DHopf = full_operator(DHopf_op, DHopf_tail_lin,DHopf_tail_nonlin, coef_nonlin);

