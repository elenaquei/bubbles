function DDH = second_der_Hopf_fullop(Xs,V)
% function DDH = second_der_Hopf_fullop(Xs,V)

ys = Xs.vector1D;
zs = squeeze(Xs.vector2D);
l1s = Xs.parameters(2);
l2s = Xs.parameters(1);
as = Xs.parameters(3);

y = V.vector1D;
z = squeeze(V.vector2D);
l1 = V.parameters(2);
l2 = V.parameters(1);
a = V.parameters(3);

y_in_z = embedding(y,z);
selected_indeces = Xs.node_time+1+(0:2*Xs.node_space)*(2*Xs.node_time+1);
z_in_y = z(selected_indeces);
ys_in_z = embedding(ys,z);
zs_in_y = zs(selected_indeces);

% 1,1 : (-K2_y^4*l2s+K2_y^2)*l1-K2_y^4*l1s*l2+4*1i*K2_y*conv_operator(l1*ys+l1s*y,1)
% 1,1 : K2^4 (-l1*l2s-l1s*l2) + K2^2*l1 + K2*4*1i*(l1*ys+l1s*y)

% 2,1 : 4*1i*K2_z*conv_operator(l1s*z+l1*zs,1)
% 2,1 : K2 * 4*1i*(l1s*z+l1*zs)

% 2,2 :  (-K2_z^4*l2s+K2_z^2)*l1-K2_z^4*l1s*l2+4*1i*K2_z*conv_operator(vec_33,1);
% 2,2 : K2^4*(l1*l2s+l1s*l2) + K2^2*l1
% vec_22 = l1*ys + l1s*y + 2*l1*as*zs+2*l1s*a*zs +2*l1s*as*z;
pol_tail1 = polynomial([-l2, 0, 0, 0, 0]);
pol_tail2 = polynomial([-l1*l2s-l1s*l2, 0, l1, 0, 0]);
zero_pol = polynomial(0);

cell_lin = cell(2);
cell_lin{1,1} = pol_tail1;
cell_lin{1,2} = zero_pol;
cell_lin{2,1} = zero_pol;
cell_lin{2,2} = pol_tail2;

%pol_nonlin = polynomial([1,0]);

% non linearities are completely in the center!
cell_nonlin = cell(2);
cell_nonlin{1,1} = zero_pol;
cell_nonlin{1,2} = zero_pol;
cell_nonlin{2,1} = zero_pol;
cell_nonlin{2,2} = zero_pol;

coef_nonlin = cell(2);
coef_nonlin{1,1} = 2*1i*y;
coef_nonlin{1,2} = 0;
coef_nonlin{2,1} = 2*1i*(l1s*z_in_y+l1*zs_in_y);
vec_22 = l1*ys_in_z + l1s*y_in_z + l1*as*zs+ l1s*a*zs +l1s*as*z;
coef_nonlin{2,2} = 2*1i*vec_22;

DDH_tail_lin = tail_matrix(cell_lin);
DDH_tail_nonlin = tail_matrix(cell_nonlin);

DDH_center = second_derHopf(Xs, V, 1,1);

DDH = full_operator(DDH_center, DDH_tail_lin,DDH_tail_nonlin, coef_nonlin);