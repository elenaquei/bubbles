function DDDH = third_derHopf_fullop(X, V, W)
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
DDDH_center = third_derHopf(X,V,W,1,1);

% define variables
y = X.vector1D;
z = squeeze(X.vector2D);
l1 = X.parameters(2);
l2 = X.parameters(1);
a = X.parameters(3);

y1 = V.vector1D;
z1 = squeeze(V.vector2D);
l11 = V.parameters(2);
l21 = V.parameters(1);
a1 = V.parameters(3);

y2 = W.vector1D;
z2 = squeeze(W.vector2D);
l12 = W.parameters(2);
l22 = W.parameters(1);
a2 = W.parameters(3);

selected_indeces = X.node_time+1+(0:2*X.node_space)*(2*X.node_time+1);

y_in_z = embedding(y,z);
y1_in_z = embedding(y1,z);
y2_in_z = embedding(y2,z);

z_in_y = z(selected_indeces);
z1_in_y = z1(selected_indeces);
z2_in_y = z2(selected_indeces);

% define the polynomial tail
pol_tail = polynomial([-l11* l22 - l12 * l21, 0, 0, 0, 0]);
zero_pol = polynomial(0);

cell_lin = cell(2);
cell_lin{1,1} = zero_pol;
cell_lin{1,2} = zero_pol;
cell_lin{2,1} = zero_pol;
cell_lin{2,2} = pol_tail;

pol_nonlin = polynomial([1,0]);

cell_nonlin = cell(2);
cell_nonlin{1,1} = pol_nonlin;
cell_nonlin{1,2} = zero_pol;
cell_nonlin{2,1} = zero_pol;
cell_nonlin{2,2} = pol_nonlin;

coef_nonlin = cell(2);
coef_nonlin{1,1} = 2 *1i*l12* z1_in_y + 2 * 1i* l11 *z2_in_y;
coef_nonlin{1,2} = 0;
coef_nonlin{2,1} =  0;
coef_nonlin{2,2} =  2 * 1i * l11 * y2_in_z + 2 * 1i * l12 * y1_in_z +... 
 l12* ( 2*a1 *1i *z + 2 *a *1i* z1) + ...
 a2* (2 *1i * l11 * z + 2* 1i *l1 * z1) + ...
 (2 *a1 *1i* l1 + 2 *a* 1i* l11)* z2;

DDDH_tail_lin = tail_matrix(cell_lin);
DDDH_tail_nonlin = tail_matrix(cell_nonlin);

% put together the full operator
DDDH = full_operator(DDDH_center, DDDH_tail_lin,DDDH_tail_nonlin, coef_nonlin);