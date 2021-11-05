% unit tests for scalar_eq with non_computable elements
% needs to be associated with the file "non_computable_func_unit_test.m"
global use_intlab
global nu
nu = 1.1;
use_intlab = 0;
Rmax = 10^-2;

load('./Kevins_code/point_candidate.mat')

% TODO: extract info from point_candidate
xi = Xi_vector([],[]);

% TODO: create new non_computable equations
n_non_comp_eqs = 2;
f = @(x) noncomputable_eqs(x, 2, [1,2,3,4,5;2,3,4,5,6]);

% not blowed up: \dot y + y - R0 h(Delay(y, tau)) * y * (1 - y)

% with polynomial embedding and blow up, with y = x + az
% F1 = - z0 + mu * x *( 1-x) * Delay( z1, tau) + mu * Delay(z2, tau) * ( -a * z0^2 + (1-2x)*z0)
% \dot z0 = F1
% \dot z1 = - p * z2 * F1 + eta1 
% \dot z2 = - a p z2 * F1 + eta2 
% scalar unknowns: psi, x, a, mu, p, eta1, eta2
% we know mu = R0 * exp(-p x)

tau = 1;

minus_F1 = '+ z0 - mu x Delay(z1,tau) + mu x^2 Delay(z1,tau) + a mu Delay(z2,tau) z0^2 - mu Delay(z2,tau) z0  + 2 x mu Delay(z2,tau) z0';
minus_F2 = '- p   z2 z0 + p   mu x Delay(z1,tau) z2 - p   mu x^2 Delay(z1,tau) z2 - p   a mu Delay(z2,tau) z0^2 z2 + p   z2 mu Delay(z2,tau) z0  - 2 p   x mu Delay(z2,tau) z0 z2 + eta1';
minus_F3 = '- a p z2 z0 + a p mu x Delay(z1,tau) z2 - a p mu x^2 Delay(z1,tau) z2 - a p a mu Delay(z2,tau) z0^2 z2 + a p z2 mu Delay(z2,tau) z0  - 2 a p x mu Delay(z2,tau) z0 z2 + eta2';

string_vector_field = strcat('psi dot z0', minus_F1, '\n psi dot z1', minus_F2, '\n psi dot z2', minus_F3);

string_vector_field = strrep(string_vector_field, 'psi' , 'l1');
string_vector_field = strrep(string_vector_field, 'x' , 'l2');
string_vector_field = strrep(string_vector_field, 'a' , 'l3');
string_vector_field = strrep(string_vector_field, 'mu' , 'l4');
string_vector_field = strrep(string_vector_field, 'p' , 'l5');
string_vector_field = strrep(string_vector_field, 'eta1' , 'l6');
string_vector_field = strrep(string_vector_field, 'eta2' , 'l7');
string_vector_field = strrep(string_vector_field, 'z0' , 'x1');
string_vector_field = strrep(string_vector_field, 'z1' , 'x2');
string_vector_field = strrep(string_vector_field, 'z2' , 'x3');
string_vector_field = strrep(string_vector_field, 'tau' , num2str(tau)); % OR - tau??

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
