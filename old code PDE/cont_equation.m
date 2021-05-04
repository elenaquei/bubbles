function [Es,const_Es] = cont_equation(x_PDE, G)
%
% INPUT 
% x_PDE     PDE_vector with 3 unknown parameters, one 1D vector and one
%           2D vector
% Es        PDE_vector or vector compatible with the previous PDE_vector,
%           continuation equation
% G         matrix of 2 lines, each line compatible with the
%           PDE_vector
% OUTPUT
% Es        elements for the continuation equation

if x_PDE.size_real ~= 3
    error('Input does not respect dimensional contraints - parameters')
end
if x_PDE.size_vector ~= 1
    error('Input does not respect dimensional contraints - vector')
end

DHopf= derivative_Hopf(x_PDE,G);

%DHopf = fit_for_inverse(DHopf_full, x_PDE.node_space);

Es = kernel(DHopf);
Es = 0 * Es; % NEW CONT EQ.
Es(3) = 1;
const_Es = - x_PDE.parameters(3);%dot( long_pyz, continuation_equation);
%const_Es = - PDE_vec2vec(x_PDE)*Es;

if issparse(Es)
    Es = full(Es);
end
%Es = adding_zero_indeces(Es,x_PDE);


% Es = kernel(DHopf);
% 
% if issparse(Es)
%     Es = full(Es);
% end
% %Es = adding_zero_indeces(Es,x_PDE);
% 
% const_Es = - PDE_vec2vec(x_PDE)*Es;