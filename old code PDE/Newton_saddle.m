 function x_zero = Newton_saddle(x_PDE, G, const_G, Es, const_Es, EDelta, const_EsDelta, tol, max_iter)
%  function x_zero = Newton_saddle(x_PDE, G, const_G, Es, const_Es, EDelta, const_EsDelta, tol, max_iter)
%
% INPUT 
% x_PDE     PDE_vector with 3 unknown parameters, one 1D vector and one
%           2D vector
%  G, const_G, Es, const_Es, EDelta, const_EsDelta        vectors compatible with the previous PDE_vector,
%           continuation equation
% tol       DEFAULT 10^-6
% max_iter  DEFAULT 30
% OUTPUT
% x_zero    solution to saddle KS via Newton method

if nargin<7 || isempty(max_iter)
    max_iter = 30;
end
if nargin<6 || isempty(tol)
    tol = 10^-6;
end


x = PDE_vec2vec(x_PDE);
y_PDE = evaluate_saddle(x_PDE, G, const_G, Es, const_Es, EDelta, const_EsDelta);
tries=0;
while norm_PDE(y_PDE)>tol && tries<max_iter
    tries=tries+1;
    y = PDE_vec2vec(y_PDE);
    derivative = derivative_saddle(x_PDE,G,Es);
    x = x.' - derivative \ y.';
    x=x.'; 
    x_PDE = PDE_vector(x, x_PDE.size_real, x_PDE.size_vector, x_PDE.node_space);
    y_PDE = evaluate_saddle (x_PDE, G, const_G,Es,const_Es, EDelta, const_EsDelta);
end
if norm_PDE(y_PDE)>tol
    error('Newton did not converge')
end
x_zero = PDE_vector(x, x_PDE.size_real, x_PDE.size_vector, x_PDE.node_space);