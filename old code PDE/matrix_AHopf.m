function A = matrix_AHopf(x_PDE, G, Es, big)
% function A = matrix_AHopf(x_PDE, G, Es, big)
%
% INPUT 
% x_PDE     PDE_vector with 3 unknown parameters, one 1D vector and one
%           2D vector
% G         matrix of 2 lines, each line compatible with the
%           PDE_vector
% Es        PDE_vector or vector compatible with the previous PDE_vector,
%           continuation equation
% big       DEFAULT 0, if big, returns expanded matrix with tails up to
%           double the number of nodes 
% OUTPUT
% DHopf     matrix, derivative of Kuramoto-shivashinki Hopf in x_PDE
% 
if nargin<4 || isempty(big)
    big =0;
end
if x_PDE.size_vector~=1
    error('wrong input')
end

DHopf = derivative_Hopf(x_PDE, G, Es);

A_small = (DHopf)^-1;
if ~big
    A = A_small;
    return
end

A = big_A(A_small,x_PDE);
if size(A,1)>10^3
    A = sparse(A);
end