function alpha = fancy_scalar_condition(xi_vec, alpha, n_eq)
% function alpha = fancy_scalar_condition(x_Xivec, alpha, n_eq)
%
% create an istance of scalar_eq with the fancy scalar condition 
%     int_[0,2pi] x(t) dot x'(t) dt = 0
% or replace such equation to the already existing istance alpha of the
% scalar_eq class.
% 
% INPUT
% x_Xivec    approximate solution
% alpha      scalar_eq where to append the new scalar eq.(DEFAULT empty)
% n_eq       int, which equation to overwright (DEFAULT 1)
% OUTPUT
% alpha      istance of scalar_eq containing the linear coefficients referring to the
%            fancy scalar condition
%                  \int_[0,2pi] x(t) dot x'(t) dt = 0
%

if nargin < 3 || isempty(n_eq)
    n_eq = 1;
end

vector_g = zeros(xi_vec.size_vector , (2*xi_vec.nodes+1));
K = -xi_vec.nodes:xi_vec.nodes;
K_inv = (2*xi_vec.nodes+1):-1:1;
for j=1:xi_vec.size_vector
    vector_g(j,:) = 1i * K .* xi_vec.vector(j,K_inv);
end


n_eqs = 1;
n_eq_vec = 0;
n_scalar = xi_vec.size_scalar;
n_vector = xi_vec.size_vector;
lin_coef = cell(3,1);
lin_coef{1} = zeros(1,n_scalar);
lin_coef{2} = zeros(1,n_vector,1+2*xi_vec.nodes);
lin_coef{2}(1,:,:) = vector_g;
lin_coef{3} = 0;
if nargin<2 || isempty(alpha)
    polynomial = polynomial_coefs(n_scalar, n_vector, 0, ...
        [], [],[],[]);
    alpha = scalar_eq(n_eqs, n_eq_vec, n_scalar, n_vector, lin_coef, polynomial);
else
    alpha = change_lin_coef(alpha,lin_coef,n_eq);
end

