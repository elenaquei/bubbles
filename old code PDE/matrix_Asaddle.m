function A = matrix_Asaddle(x_PDE, G, Es, EsDelta,big)
% function A = matrix_Asaddle(x_PDE, G, Es, EsDelta, big)
%
% INPUT 
% x_PDE     PDE_vector with 3 unknown parameters, one 1D vector and one
%           2D vector
% G         matrix of 2 lines, each line compatible with the
%           PDE_vector
% Es        PDE_vector or vector compatible with the previous PDE_vector,
%           continuation equation
% EsDelta   PDE_vector or vector compatible with the previous PDE_vector,
%           continuation equation
% big       DEFAULT 0, if big, returns expanded matrix with tails up to
%           double the number of nodes 
% _____ OR ______
% function A = matrix_Asaddle(DH, xPDE, [], [], big)
%
% INPUT 
% DH        saddle derivative derivative in x_PDE
% x_PDE     PDE_vector with 3 unknown parameters, one 1D vector and one
%           2D vector
% big       DEFAULT 0, if big, returns expanded matrix with tails up to
%           double the number of nodes 
% OUTPUT
% DHopf     matrix, derivative of Kuramoto-shivashinki Hopf in x_PDE
% 
if nargin<5 || isempty(big)
    big =0;
end

if isempty(Es)
    DH = x_PDE;
    x_PDE = G;
else
    if x_PDE.size_vector~=3
        error('wrong input')
    end
    DH = derivative_saddle(x_PDE, G, Es,EsDelta);
end

A_small = DH^-1;
if ~big
    A = A_small;
    return
end

A = big_Asaddle(A_small,x_PDE);
