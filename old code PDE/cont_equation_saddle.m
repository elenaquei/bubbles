function [Es,const_Es] = cont_equation_saddle(x_PDE, G,Es,EsDelta)
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

if x_PDE.size_real ~= 9
    error('Input does not respect dimensional contraints - parameters')
end
if x_PDE.size_vector ~= 3
    error('Input does not respect dimensional contraints - vector')
end
Es_short = Es(1:length(Es)/3);
EsDelta_short = EsDelta(1:length(EsDelta)/3);
DH= derivative_saddle (x_PDE,G,Es_short,EsDelta_short);

Es = kernel(overop2mat(DH));
const_Es = PDE_vec2vec(x_PDE)*Es;
