function x_zero = Newton_HopfKS(x_PDE, G, const_G, Es, const_Es, tol, max_iter)
% function x_zero = Newton_HopfKS(x_PDE, G, Es)
%
% INPUT 
% x_PDE     PDE_vector with 3 unknown parameters, one 1D vector and one
%           2D vector
% Es        PDE_vector or vector compatible with the previous PDE_vector,
%           continuation equation
% G         matrix of 2 lines, each line compatible with the
%           PDE_vector
% Es        (DEFAULT empty) elements for the continuation equation
% OUTPUT
% x_zero    solution to Hopf KS via Newton method
global display

if nargin<7 || isempty(max_iter)
    max_iter = 30;
end
if nargin<6 || isempty(tol)
    tol = 10^-9;
end


x = PDE_vec2vec(x_PDE);
y_PDE = evaluate_HopfKS(x_PDE, G, const_G, Es, const_Es);
tries=0;
Remainder = zeros(max_iter);
while norm_PDE(y_PDE)>tol && tries<max_iter
    Remainder(tries+1) = norm_PDE(y_PDE);
    tries=tries+1;
    y = PDE_vec2vec(y_PDE);
    derivative = derivative_Hopf(x_PDE,G,Es,[],0);
    %DH = Hopf_derivative_FD(x_PDE, G, Es);
    %inv_der = adding_zero_indeces(fit_for_inverse(derivative,x_PDE)^-1,x_PDE);
    %rcond(DH);
    %if any(abs(derivative-DH)>10^-4)
    %    error('wrong der')
    %end
    %norm(y)
    x = x.' -  derivative\ y.';
    x=x.'; 
    x_PDE = PDE_vector(x, x_PDE.size_real, x_PDE.size_vector, x_PDE.node_space);
    y_PDE = evaluate_HopfKS(x_PDE, G, const_G,Es,const_Es);
end
if norm_PDE(y_PDE)>tol
    Remainder(tries+1) = norm_PDE(y_PDE);
    error('Newton did not converge, remainder %1.4e',Remainder(tries+1))
end
Remainder = Remainder(1:tries+1);
if display >1
    fprintf('Netwon remainder: %1.4d\n',Remainder(end))
end

x_zero = PDE_vector(x, x_PDE.size_real, x_PDE.size_vector, x_PDE.node_space);