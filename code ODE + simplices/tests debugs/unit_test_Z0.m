% unit tests for scalar_eq with non_computable elements
% needs to be associated with the file "non_computable_func_unit_test.m"
global use_intlab
global nu
nu = 1.1;
use_intlab = 0;
Rmax = 10^-2;
x1 = [0,-1i,1,+1i,0];
x2 = [1,2,3,2,1] + [0,1,0,-1,0]*1i;
x3 = [2,-1,0,-1,2] + [1,1,0,-1,-1]*1i;



l1 = -2;
l2 = 1.2;
h = 10^-3;
xi_0 = Xi_vector([l1,l2], [x1;x2;x3]);
xi_1 = Xi_vector([l1,l2]+ rand(1,2)*h, [x1;x2;x3]+ rand(3,length(x1))*h);
xi_2 = Xi_vector([l1,l2]+ rand(1,2)*h, [x1;x2;x3]+ rand(3,length(x1))*h);

n_non_comp_eqs = 2;
f = @(x) noncomputable_eqs(x, 2, [1,2,3,4,5;2,3,4,5,6]);

string_vector_field = 'dot x1 + l1 * Delay(x1,3.14) * Delay(x2,4) + x2^2 * x3 \n dot x2 + l2*Delay(x2,3) * Delay(x1,2) *x3 - x1 \n dot x3 + x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,2,3);

n_scal = xi_0.size_scalar;
n_vec = xi_0.size_vector;

polynomial = polynomial_coefs(n_scal, n_vec, 0, ...
    [], [],[],[]);
scalar_equation = scalar_eq(0, 0, n_scal,...
    n_vec, cell(3,1), polynomial, n_non_comp_eqs, f);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

use_intlab = 1;
Z0 = Z0_delay_simplex(full_zero_finding_problem, xi_0, xi_1, xi_2);

Y = Y_delay_simplex(full_zero_finding_problem, xi_0, xi_1, xi_2);

Z1 = Z1_delay_simplex(full_zero_finding_problem, xi_0, xi_1, xi_2);

Z2 = Z2_delay_simplex(full_zero_finding_problem, xi_0, xi_1, xi_2, Rmax);
use_intlab = 0;

function [F, dxF, dxF_mat] = noncomputable_eqs(x, n_eqs, constants)
% function noncomputable_eqs(x, n_eqs)
global use_intlab
if nargin > 2
    if size(constants,1)~=n_eqs
        constants = constants.';
    end
    if size(constants,1)~=n_eqs || size(constants,2)~= x.size_vector + x.size_scalar
        error('constants has the wrong size');
    end
end
if nargin < 2
    n_eqs = 1;
end
if nargin < 3
    constants = ones(n_eqs, x.size_vector + x.size_scalar);
end

sum_x = zeros(x.size_vector + x.size_scalar,1);
if isintval(x.scalar) || isintval(x.vector) || use_intlab
    sum_x = intval(sum_x);
end
for j = 1: x.size_scalar
    sum_x(j) = x.scalar(j);
end
for j = 1:x.size_vector
    sum_x(x.size_scalar + j) = sum(x.vector(j,:));
end

F = constants * sum_x;
if nargout == 1
    return
end
dxF = constants;
if nargout == 2
    return
end
x_one_vec = 1 + 0*x.vector(1,:);

dxF_mat = [constants(:,1:x.size_scalar), kron(constants(:,(x.size_scalar+1):end), x_one_vec)];
return
end
