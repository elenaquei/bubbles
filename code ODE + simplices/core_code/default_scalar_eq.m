function alpha = default_scalar_eq(xi_vec,n_component, val, alpha, n_eq)
% function alpha = default_scalar_eq(xi_vec,n,val)
%
% creates a default scalar equation or append it to the already existing 
% scalar_eq alpha. It adds just a linear equation summing
% on the n_component-th component of xi_vec and keeping it constant or equal to
% the given value val.
%
% INPUTS
% xi_vec            Xi_vector, approximation of the wated solution
% n_component       which component to take (DEFAULT 1)
% val               which value to give, must be real
% alpha             scalar_eq where to append the new scalar eq.(DEFAULT empty)
% n_eq              int, which equation to overwright (DEFAULT 1)
% OUTPUTS
% alpha             scalar_eq with the appendind equation


if nargin<2 || isempty(n_component)
    n_component=1;
end
if nargin<3 || isempty(val)
    val = real(sum(xi_vec.vector(n_component,:)));
end
if nargin<5 || isempty(n_eq)
    n_eq = 1;
end

n_eqs = 1;
n_eq_vec = 0;
n_scalar = xi_vec.size_scalar;
n_vector = xi_vec.size_vector;
lin_coef = cell(3,1);
lin_coef{1} = zeros(1,n_scalar);
lin_coef{2} = zeros(1,n_vector,1+2*xi_vec.nodes);
lin_coef{2}(1,n_component,:)=1;
lin_coef{3} = - val;

if nargin < 4 || isempty(alpha)
    polynomial = polynomial_coefs(n_scalar, n_vector, 0, ...
        [], [],[],[]);
    alpha = scalar_eq(n_eqs, n_eq_vec, n_scalar, n_vector, lin_coef, polynomial);
else
    alpha = change_lin_coef(alpha,lin_coef,n_eq);
end

end
%end DEFAULT_SCALAR_EQ