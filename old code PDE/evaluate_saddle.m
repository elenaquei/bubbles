function y_PDE = evaluate_saddle(xlong_PDE, G, const_G,Es,const_Es,EDelta,const_EsDelta, flagFullSize)
% function y_PDE = evaluate_saddle(x_PDE, G, const_G,Es,const_Es, flagFullSize)
%
% INPUT
% x_PDE         PDE_vector with
%                   size_real 9
%                   size_vector 3
% G             compatible with a PDE_vector of size_real 3 and size_vector
%                    1
% const_G, Es, const_Es as G
% flagFullSize  bool, DEFAULT 0
% OUTPUT
% y_PDE         PDE_vector with
%                   size_real 9
%                   size_vector 3
% if flagFullSize y_PDE has double the nodes as x_PDE, otherwise they are
% of hte same lenght

% if nargin>7 && flagFullSize
%     xlong_PDE = reshape(xlong_PDE, xlong_PDE.node_space*2, xlong_PDE.node_time*2);
% end
if nargin==7
    flagFullSize=0;
end
if size(G,1)~=4
    error('Size G incompatible')
end
if xlong_PDE.size_vector~=3 || xlong_PDE.size_real~=9
    error('x_PDE incompatible')
end

[x_PDE, v_PDE, w_PDE] = split(xlong_PDE);
v_vec_small = v_PDE;
if flagFullSize
    v_PDE = reshape(v_PDE, v_PDE.node_space*2, v_PDE.node_time*2);
    w_PDE = reshape(w_PDE, w_PDE.node_space*2, w_PDE.node_time*2);
end

x_vec = PDE_vec2vec(x_PDE);
v_vec = PDE_vec2vec(v_PDE);
%w_vec = PDE_vec2vec(w_PDE);

DH = derivative_Hopf(x_PDE,G,Es,flagFullSize);
y1_PDE = evaluate_HopfKS(x_PDE, G, const_G,Es,const_Es, flagFullSize);

y2_PDE = PDE_vector(DH*PDE_vec2vec(v_PDE).',v_PDE);
y2_PDE.parameters(1)=y2_PDE.parameters(1) + compute_linear_func(x_vec, EDelta(:).',const_EsDelta);

if flagFullSize
    x_PDE = reshape(x_PDE, x_PDE.node_space*2, x_PDE.node_time*2);
end

%y2_PDE = PDE_vector(y2_vec, 3, 1, x_PDE.node_space);


y3_PDE = second_derHopf(x_PDE,v_PDE,0)*v_PDE+...
    DH*w_PDE;
y3_PDE.parameters(1) = y3_PDE.parameters(1) + ...
    2 * compute_linear_func(v_vec_small, EDelta(:).',0*const_EsDelta);

%y3_PDE = PDE_vector(y3_vec, 3, 1, x_PDE.node_space);

y_PDE = merge(y1_PDE, y2_PDE, y3_PDE);








