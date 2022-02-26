function DH_nonlin = derivative_Hopf_fullop_minus_Adag(X)


G = zeros(4,length(X));
Es = zeros(length(X),1);
DHopf_op = derivative_Hopf(X, G, Es, 1,1);
DHopf_op = tail(DHopf_op,X.node_space,X.node_time);
DHopf_op = non_diag(DHopf_op);

lambda = X.parameters(2);
% nu = X.parameters(1);
a = X.parameters(3);

y = X.vector1D;
z = squeeze(X.vector2D);

y_in_z = embedding(y,z);
selected_indeces = X.node_time+1+(0:2*X.node_space)*(2*X.node_time+1);
z_in_y = z(selected_indeces);

zero_pol = polynomial(0);

cell_lin = cell(2);
cell_lin{1,1} = zero_pol;
cell_lin{1,2} = zero_pol;
cell_lin{2,1} = zero_pol;
cell_lin{2,2} = zero_pol;

pol_nonlin = polynomial([1,0]);

cell_nonlin = cell(2);
cell_nonlin{1,1} = pol_nonlin;
cell_nonlin{1,2} = zero_pol;
cell_nonlin{2,1} = pol_nonlin;
cell_nonlin{2,2} = pol_nonlin;

coef_nonlin = cell(2);
coef_nonlin{1,1} = 2*1i*lambda*y;
coef_nonlin{1,2} = 0;
coef_nonlin{2,1} = 2*lambda*1i*z_in_y;
coef_nonlin{2,2} = 2*lambda*1i*(y_in_z + a*z);

DHopf_tail_lin = tail_matrix(cell_lin);
DHopf_tail_nonlin = tail_matrix(cell_nonlin);

DH_nonlin = full_operator(DHopf_op, DHopf_tail_lin,DHopf_tail_nonlin, coef_nonlin);


