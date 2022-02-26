function [beta,x_dot_1] = continuation_eq(alpha,x1,old_x_dot,step_size)
% function beta = continuation_eq(alpha,x1,old_x_dot,step_size)
%
% INPUT
% alpha         full_problem, not square
% x1            Xi_vector, solution of alpha
% old_x_dot     (DEFAULT = []) dot derivative of previous step
% step_size     (DEFAULT = []) for check that derivative coherent
% OUTPUT
% beta          full_problem, square with continuation equation

beta = alpha;

% compute kernel and set it real
x_dot_1 = kernel(derivative_to_matrix(derivative(alpha, x1,0)));
[~,index] = max(abs(real(x_dot_1(1:x1.size_scalar))));
angle = atan( imag(x_dot_1(index))/real(x_dot_1(index)));
x_dot_1 = exp( - 1i * angle) * x_dot_1;
x_dot_1 = Xi_vec2vec(symmetrise(vec2Xi_vec(x_dot_1,x1)));
if nargin>= 3  && dot(x_dot_1,old_x_dot)<0
    x_dot_1=-x_dot_1;
end
const = -Xi_vec2vec(x1).'*conj(x_dot_1);

% addition of continuation equation
beta.scalar_equations = change_lin_coef_vector(alpha.scalar_equations,[conj(x_dot_1);const],alpha.scalar_equations.number_equations_lin+1);


if nargin>3 && norm(old_x_dot-x_dot_1) > 100*step_size
    warning('Norm derivatives in x0 and x1 quite big')
end