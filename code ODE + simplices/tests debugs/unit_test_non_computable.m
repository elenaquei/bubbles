function unit_test_non_computable()
% unit tests for scalar_eq with non_computable elements
% needs to be associated with the file "non_computable_func_unit_test.m"
x1 = [0,-1i,1,+1i,0];
x2 = [1,2,3,2,1] + [0,1,0,-1,0]*1i;
x3 = [2,-1,0,-1,2] + [1,1,0,-1,-1]*1i;

l1 = -2;
l2 = 1.2;
x = Xi_vector([l1,l2], [x1;x2;x3]);
one_vec = [0,0,1,0,0];
mat_K = diag(-2:2);
mat_bigK = diag(-4:4);
mat_verybigK = diag(-6:6);
conv_mat = @(x) conj(toeplitz(x((length(x)+1)/2:end)));

conv_mat_big =@(x) conj(toeplitz([x((length(x)+1)/2:end),0*x((length(x)+3)/2:end)]));
Delay = @(x, tau) exp(1i*l1*(-(length(x)-1)/2:(length(x)-1)/2)*tau).*x;

n_non_comp_eqs = 2;
f = @(x) noncomputable_eqs(x, 2, [1,2,3,4,5;2,3,4,5,6]);
% [a,b,c]= f(x);

string_vector_field = 'dot x1 + l1 * Delay(x1,3.14) * Delay(x2,4) + x2^2 * x3 \n dot x2 + l2*Delay(x2,3) * Delay(x1,2) *x3 - x1 \n dot x3 + x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,2,3);
n_scal = x.size_scalar;
n_vec = x.size_vector;

polynomial = polynomial_coefs(n_scal, n_vec, 0, ...
    [], [],[],[]);
scalar_equation = scalar_eq(0, 0, n_scal,...
    n_vec, cell(3,1), polynomial, n_non_comp_eqs, f);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_sol = apply(full_zero_finding_problem, x);
[hand_sol, hand_small_der, hand_big_der] = f(reshape(x, x.nodes*3));

auto_der = derivative(full_zero_finding_problem, x);

% string_vector_field = 'dot x1 + l1 * Delay(x1,3.14) * Delay(x2,4) + x2^2 * x3 \n dot x2 + l2*Delay(x2,3) * Delay(x1,2) *x3 - x1 \n dot x3 + x1';
auto_mat_der = derivative_to_matrix(auto_der);
one_vec_big = conv(one_vec,one_vec);

dot_der = kron([1,0,0; 0,1,0;0,0,1], mat_verybigK * 1i);
mat_conv = [ 0*conv_mat_big(conv(Delay(x1,3.14),one_vec_big)), conv_mat_big(conv(2*conv(x2,x3),one_vec)) , conv_mat_big(conv(conv(x2,x2),one_vec))
    conv_mat_big(conv(-one_vec_big,one_vec)), conv_mat_big(conv(0*one_vec_big,one_vec)) , l2*conv_mat_big(conv(conv(Delay(x2,3),Delay(x1,2)),one_vec))
    conv_mat_big(conv(one_vec_big,one_vec)) , conv_mat_big(conv(0*one_vec_big,one_vec)) , conv_mat_big(conv(one_vec_big*0,one_vec)) ];
d_lambda = [conv(one_vec,conv(Delay(x1,3.14), Delay(x2,4))), 0*conv(one_vec_big, x2), conv(one_vec_big, 0*one_vec);
    0*conv(one_vec,conv(Delay(x2,4), Delay(x1,3.14))), conv(x3, conv(Delay(x2,3), Delay(x1,2))), conv(one_vec_big, 0*one_vec)].';

delay_mat11 = conv_mat_big(l1*conv(Delay(x2,4),one_vec_big)) *  expm(1i*3.14*l1 * mat_verybigK);
delay_mat12 = conv_mat_big(l1*conv(Delay(x1,3.14),one_vec_big)) *  expm(1i*4*l1 * mat_verybigK);
delay_mat21 = conv_mat_big(l2*conv(conv(Delay(x2,3),x3),one_vec)) *  expm(1i*2*l1 * mat_verybigK);
delay_mat22 = conv_mat_big(l2*conv(Delay(x1,2),conv(x3,one_vec))) *  expm(1i*3*l1 * mat_verybigK);
delay_mat = [delay_mat11, delay_mat12, 0*delay_mat11
    delay_mat21, delay_mat22, 0*delay_mat22
    0*delay_mat22, 0*delay_mat22, 0*delay_mat22];

hand_solution = [hand_big_der;
    d_lambda, dot_der + delay_mat + mat_conv];

max_error = norm(hand_solution - auto_mat_der);

if max(max_error) > 10^-10
    error('Test 4 derivative delay - 3D to 3D')
end

end

function [F, dxF, dxF_mat] = noncomputable_eqs(x, n_eqs, constants)
% function noncomputable_eqs(x, n_eqs)

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
for j = 1: x.size_scalar
    sum_x(j) = x.scalar(j);
end
for j = 1:x.size_vector
    sum_x(x.size_scalar + j) = sum(x.vector(j,:));
end

F = constants * sum_x;
dxF = constants;

x_one_vec = 1 + 0*x.vector(1,:);

dxF_mat = [constants(:,1:x.size_scalar), kron(constants(:,(x.size_scalar+1):end), x_one_vec)];
return
end
