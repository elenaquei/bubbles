function DDDDH = fourth_derHopf_fullop(X, V, W, Z)
% INPUT
% X     PDE_vector
% V     PDE_vector, first direction
% w     PDE_vector, second direction
% flag  bool
% OUTPUT
% DDDDHxDxD  full_operator

if nargin == 2
    W = V;
end

% construct the overoperator in the center
DDDH_center = fourth_derHopf(X,V,W,Z,1);

y = X.vector1D;
z = squeeze(X.vector2D);
l1 = X.parameters(1);
l2 = X.parameters(2);
a = X.parameters(3);

y1 = V.vector1D;
z1 = squeeze(V.vector2D);
l11 = V.parameters(1);
l21 = V.parameters(2);
a1 = V.parameters(3);

y2 = V.vector1D;
z2 = squeeze(V.vector2D);
l12 = V.parameters(1);
l22 = V.parameters(2);
a2 = V.parameters(3);

y3 = V.vector1D;
z3 = squeeze(V.vector2D);
l13 = V.parameters(1);
l23 = V.parameters(2);
a3 = V.parameters(3);

selected_indeces = X.node_time+1+(0:2*X.node_space)*(2*X.node_time+1);

y_in_z = embedding(y,z);
y1_in_z = embedding(y1,z);
y2_in_z = embedding(y2,z);
y3_in_z = embedding(y3,z);

z_in_y = z(selected_indeces);
z1_in_y = z1(selected_indeces);
z2_in_y = z2(selected_indeces);
z3_in_y = z3(selected_indeces);

% define the polynomial tail
%pol_tail = polynomial([-l11* l22 - l12 * l21, 0, 0, 0, 0]);
zero_pol = polynomial(0);

cell_lin = cell(2);
cell_lin{1,1} = zero_pol;
cell_lin{1,2} = zero_pol;
cell_lin{2,1} = zero_pol;
cell_lin{2,2} = zero_pol;

pol_nonlin = polynomial([1,0]);

cell_nonlin = cell(2);
cell_nonlin{1,1} = zero_pol;
cell_nonlin{1,2} = zero_pol;
cell_nonlin{2,1} = zero_pol;
cell_nonlin{2,2} = pol_nonlin;

coef_nonlin = cell(2);
coef_nonlin{1,1} = 0;
coef_nonlin{1,2} = 0;
coef_nonlin{2,1} =  0;
coef_nonlin{2,2} =  l13 *(2 *a2* 1i *z1 + 2* a1 *1i *z2) +... 
 a3* (2* 1i *l12 *z1 + 2 *1i *l11 *z2) + ...
 (2* a2* 1i *l11 + 2 *a1* 1i *l12)* z3;

DDDDH_tail_lin = tail_matrix(cell_lin);
DDDDH_tail_nonlin = tail_matrix(cell_nonlin);

% put together the full operator
DDDDH = full_operator(DDDH_center, DDDDH_tail_lin,DDDDH_tail_nonlin, coef_nonlin);