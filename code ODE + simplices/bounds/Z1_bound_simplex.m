function [Z1vector,new_Z1,Z1_s]=Z1_bound_simplex(A0,A1,A2,xBar0,xBar1,xBar2,alpha0,alpha1,alpha2,Adagger_delta1,Adagger_delta2,old_Z1_x0,old_Z1_x1)
% function [Z1vector,new_Z1,Z1_s]=Z1_bound_simplex(A0,A1,A2,xBar0,xBar1,xBar2,alpha0,alpha1,alpha2,Adagger_delta1,Adagger_delta2,old_Z1_x0,old_Z1_x1)
% 
% if given, old_Z1 needs to have the following structure:
%   old_Z1.mat (matrix)
%   old_Z1.nodes_col (positive integer)
%   old_Z1.nodes_row (positive integer)
%
% Z1_s        Z1+ s^2 Z1Delta, dependency on the stepsize of Z1 (rough bound),
%             saved as standard vector of three elements [0,Z1Delta,0, Z1]


global use_intlab
global nu
global talkative

N=xBar0.size_vector;
M=xBar0.size_scalar;

if nargin>11 && ~isempty(old_Z1_x0)
    Z_mat0=old_Z1.mat;
    nodes_col0=old_Z1.nodes_col;
    nodes_row0=old_Z1.nodes_row;
else
    [~,Z_mat0,nodes_col0,nodes_row0]=Z1_bound(A0,xBar0,alpha0);
end
if nargin>12 && ~isempty(old_Z1_x1)
    Z_mat1=old_Z1.mat;
else
    [~,Z_mat1]=Z1_bound(A1,xBar1,alpha1);
end

[~,Z_mat2,nodes_col2,nodes_row2]=Z1_bound(A2,xBar2,alpha2);
new_Z1.mat=Z_mat2;
new_Z1.nodes_col=nodes_col2;
new_Z1.nodes_row=nodes_row2;

if nodes_col2~=nodes_col0 || nodes_row2~=nodes_row0
    [~,Z_mat0,nodes_col0,nodes_row0]=Z1_bound(A0,xBar0,alpha0);
    if nodes_col1~=nodes_col0 || nodes_row1~=nodes_row0
        error('not matching dimensions');
    end
end

if use_intlab
    Z1_extrema=norm_Ximat(max(max(Z_mat0.sup,Z_mat1.sup),Z_mat2.sup),M,N,nodes_col0,nodes_row0);
else
    Z1_extrema=norm_Ximat(max(max(Z_mat0,Z_mat1),Z_mat2),M,N,nodes_col0,nodes_row0);
end


xDelta1=norm_Xi_vector(xBar0-xBar1,nu)/2;
xDelta2=norm_Xi_vector(xBar0-xBar2,nu)/2;
xDelta = max(xDelta1, xDelta2);
xBarcenter = (xBar2 + xBar1 + xBar0) ./ 3;

DDDP = Function_third_derivative(alpha0,xBarcenter,xDelta);
DDP=Function_second_derivative(alpha0,xBarcenter,xDelta);

% addition for the new continuation equation and linear in s scalar eqs
for i = 1:alpha1.scalar_equations.number_equations_lin
    coef_1 = extract_lin_coef(alpha1.scalar_equations, i);
    coef_0 =  extract_lin_coef(alpha0.scalar_equations, i);
    coef_Delta1 = coef_1-coef_0;
    DDP(i) = max(norm(vec2Xi_vec(coef_Delta1,xBar0)));
end

%
% term : A_D * partial_s H_s(x_s) x_D
% 
partial_sH_sx_D1 = zeros(length(DDP),1);
partial_sH_sx_D2 = zeros(length(DDP),1);
vec_Delta1 = Xi_vec2vec(xBar0-xBar1);
vec_Delta2 = Xi_vec2vec(xBar0-xBar2);
if isintval(vec_Delta1)
    partial_sH_sx_D1 = intval(partial_sH_sx_D1);
    partial_sH_sx_D2 = intval(partial_sH_sx_D2);
end
for i = 1:alpha1.scalar_equations.number_equations_lin
    coef_2 = extract_lin_coef(alpha2.scalar_equations, i);
    coef_1 = extract_lin_coef(alpha1.scalar_equations, i);
    coef_0 =  extract_lin_coef(alpha0.scalar_equations, i);
    const2 = alpha2.scalar_equations.linear_coef{3}(i);
    const1 = alpha1.scalar_equations.linear_coef{3}(i);
    const0 =  alpha0.scalar_equations.linear_coef{3}(i);
    coef_Delta1 = coef_1-coef_0;
    const_Delta1 = const1-const0;
    partial_sH_sx_D1(i) = abs(dot(coef_Delta1,vec_Delta1)+const_Delta1);
    coef_Delta2 = coef_2-coef_0;
    const_Delta2 = const2-const0;
    partial_sH_sx_D2(i) = abs(dot(coef_Delta2,vec_Delta2)+const_Delta2);
end
if isintval(vec_Delta1)
    partial_sH_sx_D1 = sup(partial_sH_sx_D1);
    partial_sH_sx_D2 = sup(partial_sH_sx_D2);
end

% creating interval matrix for As
if isintval(A0)
    As = infsup (   ...
    min(min(inf(real(A0)),inf(real(A1))),inf(real(A2))), ...
    max( max(sup(real(A0)),sup(imag(A1))), sup(real(A2))))+...
        1i*infsup (   ...
    min(min(inf(imag(A0)),inf(imag(A1))),inf(imag(A2))), ...
    max( max(sup(imag(A0)),sup(imag(A1))), sup(imag(A2))));
else
    As = infsup( min(min(real(A0),real(A1)), real(A2)),max(max(real(A0),real(A1)),real(A2)) )+...
        1i*infsup( min(min(imag(A0),imag(A1)), imag(A2)),max(max(imag(A0),imag(A1)),imag(A2)));
end

As_ij= component_matrix_norm(As,xBar0);
Adiff1_ij= component_matrix_norm(A1-A0,xBar0); % A_\Delta
Adiff2_ij= component_matrix_norm(A2-A0,xBar0); % A_\Delta

sum_of_derivatives_s1s1=As_ij*DDDP*xDelta1*xDelta1+2*Adiff1_ij*DDP*xDelta1 + 2*Adiff1_ij*partial_sH_sx_D1;
sum_of_derivatives_s1s2=As_ij*DDDP*xDelta1*xDelta2+2*Adiff1_ij*DDP*xDelta2 + 2*Adiff1_ij*partial_sH_sx_D2;
sum_of_derivatives_s2s1=As_ij*DDDP*xDelta2*xDelta1+2*Adiff2_ij*DDP*xDelta1 + 2*Adiff2_ij*partial_sH_sx_D1;
sum_of_derivatives_s2s2=As_ij*DDDP*xDelta2*xDelta2+2*Adiff2_ij*DDP*xDelta2 + 2*Adiff2_ij*partial_sH_sx_D2;

sum_of_derivatives = max( max(sum_of_derivatives_s1s1,sum_of_derivatives_s1s2 ), ...
    max(sum_of_derivatives_s2s1,sum_of_derivatives_s2s2));
%if use_intlab
    OnesNM=(ones(N+M,1));
%else
%    OnesNM=ones(N+M,1);
%end
%addition_to_max=1/8*(2*component_matrix_norm((A1-A0)*Adagger_delta,xBar1)*OnesNM+sum_of_derivatives);
Z1DDelta=sum_of_derivatives/8;
Z1Delta_s1s1=1/8*(2*component_matrix_norm((A1-A0)*Adagger_delta1,xBar1)*OnesNM);
Z1Delta_s1s2=1/8*(2*component_matrix_norm((A1-A0)*Adagger_delta2,xBar1)*OnesNM);
Z1Delta_s2s1=1/8*(2*component_matrix_norm((A2-A0)*Adagger_delta1,xBar1)*OnesNM);
Z1Delta_s2s2=1/8*(2*component_matrix_norm((A2-A0)*Adagger_delta2,xBar1)*OnesNM);
Z1Delta = max( max(Z1Delta_s1s1,Z1Delta_s1s2 ), ...
    max(Z1Delta_s2s1,Z1Delta_s2s2));
%addition_to_max=Z1Delta+Z1DDelta;


Z1_s= [0*Z1_extrema,Z1DDelta+Z1Delta,0*Z1_extrema, Z1_extrema];
Z1vector = sum(Z1_s,2);

% error handeling 
if any(Z1vector>1)
    fprintf('Z1_simplex computed,\n');fprintf(' %d\n',max(Z1vector));
    %error('Z1_cont is bigger than 1, no interval found')
elseif talkative>1
    fprintf('Z1_simplex computed, %d\n',max(Z1vector));
elseif talkative>0
    fprintf('Z1_simplex computed, %d\n',max(Z1vector));
end

