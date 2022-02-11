% unit test delay functionalities
% PARTS:
% - test 1D vector field (linear, polynomial, delay linear, delay polynomial,
% sums)
% - test higher D vector field
% - test derivatives


%% PART 1: 1D vector field in 1 variable

x1 = [0,-1i,1,+1i,0];
l1 = 1;
x = Xi_vector(l1, x1);
string_vector_field = 'dot x1 + l1 * x1^2 ';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,1);

auto_sol = apply(vector_field, x);
auto_sol = verifyfft_in(auto_sol, 1).';
% STILL RETURNS THE IFFT - but I don't seem to be able to update this... whyyyyy??

K = (-2:2);
exact_sol = 1i*K.*x1 + l1 * conv(x1,x1,'same'); %+ 2*x1

max_error = max(abs(auto_sol(7:11) - exact_sol));
if max_error > 10^-10
    error('Test 1 - 1D')
end


x1 = [0,-1i,1,+1i,0];
l1 = 1;
x = Xi_vector(l1, x1);
string_vector_field = 'dot x1 + l1 * x1^2 + 2 * x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,1);

auto_sol = apply(vector_field, x);
auto_sol = verifyfft_in(auto_sol, 1).';
% STILL RETURNS THE IFFT - but I don't seem to be able to update this... whyyyyy??

K = (-2:2);
exact_sol = 1i*K.*x1 + l1 * conv(x1,x1,'same') + 2*x1;

max_error = max(abs(auto_sol(7:11) - exact_sol));
if max_error > 10^-10
    error('Test 2 - 1D')
end



x1 = [0,-1i,1,+1i,0];
l1 = 1.5;
x = Xi_vector(l1, x1);
string_vector_field = 'dot x1 + l1 * Delay(x1,3.14)*x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,1);

auto_sol = apply(vector_field, x);
auto_sol = verifyfft_in(auto_sol, 1).';
% STILL RETURNS THE IFFT - but I don't seem to be able to update this... whyyyyy??

K = (-2:2);
K_big = -4:4;
x1_big = [0,0,x1,0,0];
exact_sol = 1i*K_big.*x1_big + l1 *  conv(exp(1i * l1 * K * 3.14) .* x1,x1);

max_error = max(abs(auto_sol(5:end-3) - exact_sol));
if max_error > 10^-10
    error('Test 3 - 1D')
end



x1 = [0,-1i,1,+1i,0];
l1 = 1;
x = Xi_vector(l1, x1);
string_vector_field = 'Delay(x1,3.14)*x1^2';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,1);

auto_sol = apply(vector_field, x);
auto_sol = verifyfft_in(auto_sol, 1).';
% STILL RETURNS THE IFFT - but I don't seem to be able to update this... whyyyyy??

K = (-2:2);
exact_sol = conv(conv(exp(1i * l1 * K * 3.14) .* x1,x1,'same'),x1,'same');

max_error = max(abs(auto_sol(7:11) - exact_sol));
if max_error > 10^-10
    error('Test 4 - 1D')
end


x1 = [0,-1i,1,+1i,0];
l1 = 1.5;
x = Xi_vector(l1, x1);
string_vector_field = 'dot x1 + l1 * Delay(x1,3.14)*x1^2';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,1);

auto_sol = apply(vector_field, x);
auto_sol = verifyfft_in(auto_sol, 1).';
% STILL RETURNS THE IFFT - but I don't seem to be able to update this... whyyyyy??

K = (-2:2);
exact_sol =1i*K.*x1 + l1 * conv(conv(exp(1i * l1 * K * 3.14) .* x1,x1,'same'),x1,'same');

max_error = max(abs(auto_sol(7:11) - exact_sol));
if max_error > 10^-10
    error('Test 5 - 1D')
end


x1 = [0,-1i,1,+1i,0];
l1 = 1.5;
x = Xi_vector(l1, x1);
string_vector_field = 'dot x1 + l1 * Delay(x1,3.14)*x1^2 + 2 * l1 * x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,1);

auto_sol = apply(vector_field, x);
auto_sol = verifyfft_in(auto_sol, 1).';
% STILL RETURNS THE IFFT - but I don't seem to be able to update this... whyyyyy??

K = (-2:2);
exact_sol =1i*K.*x1 + l1 * conv(conv(exp(1i * l1 * K * 3.14) .* x1,x1,'same'),x1,'same') + 2 * l1 * x1;

max_error = max(abs(auto_sol(7:11) - exact_sol));
if max_error > 10^-10
    error('Test 6 - 1D')
end



x1 = [0,-1i,1,+1i,0];
l1 = 1.5;
x = Xi_vector(l1, x1);
string_vector_field = 'dot x1 + l1 * Delay(x1,3.14)^2*x1 + 2 * l1 * x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,1);

auto_sol = apply(vector_field, x);
auto_sol = verifyfft_in(auto_sol, 1).';

K = (-2:2);
exact_sol =1i*K.*x1 + l1 * conv(conv(exp(1i * l1 * K * 3.14) .* x1,...
    exp(1i * l1 * K * 3.14) .* x1,'same'),x1,'same') + 2 * l1 * x1;

max_error = max(abs(auto_sol(7:11) - exact_sol));
if max_error > 10^-10
    error('Test 7 - 1D')
end


disp('The evaluation of a 1D vector field with delay and convolutions is successful')


%% PART 2: higher dimensions vector field
%
% subpart 1: higher dim vector, 1D vector field


x1 = [0,-1i,1,+1i,0];
x2 = [1,2,3,2,1] + [0,1,0,-1,0]*1i;
x3 = [2,-1,0,-1,2] + [1,1,0,-1,-1]*1i;

l1 = 1;
x = Xi_vector(l1, [x1;x2;x3]);
string_vector_field = 'dot x1 + l1 * x1 + x2 + x3 ';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

auto_sol = apply(vector_field, x);
auto_sol = verifyfft_in(auto_sol, 1).';
% STILL RETURNS THE IFFT - but I don't seem to be able to update this... whyyyyy??

K = (-2:2);
exact_sol = 1i*K.*x1 + l1 * x1 + x2 + x3;

max_error = max(abs(auto_sol(3:end-1) - exact_sol));
if max_error > 10^-10
    error('Test 1 - 3D to 1D')
end


string_vector_field = 'dot x1 + l1 * x1^2 + 2 * x1 * x2';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

auto_sol = apply(vector_field, x);
auto_sol = verifyfft_in(auto_sol, 1).';
% STILL RETURNS THE IFFT - but I don't seem to be able to update this... whyyyyy??

K = (-2:2);
exact_sol = 1i*K.*x1 + l1 * conv(x1,x1,'same') + 2*conv(x1,x2,'same');

max_error = max(abs(auto_sol(7:11) - exact_sol));
if max_error > 10^-10
    error('Test 2 - 3D to 1D')
end


string_vector_field = 'dot x1 + l1 * Delay(x1,3.14)*x2';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

auto_sol = apply(vector_field, x);
auto_sol = verifyfft_in(auto_sol, 1).';
% STILL RETURNS THE IFFT - but I don't seem to be able to update this... whyyyyy??

K = (-2:2);
K_big = -4:4;
x1_big = [0,0,x1,0,0];
exact_sol = 1i*K_big.*x1_big + l1 *  conv(exp(1i * l1 * K * 3.14) .* x1,x2);

max_error = max(abs(auto_sol(5:end-3) - exact_sol));
if max_error > 10^-10
    error('Test 3 - 3D to 1D')
end



string_vector_field = 'Delay(x1,3.14)*x2*x3';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

auto_sol = apply(vector_field, x);
auto_sol = verifyfft_in(auto_sol, 1).';
% STILL RETURNS THE IFFT - but I don't seem to be able to update this... whyyyyy??

K = (-2:2);
exact_sol = conv(conv(exp(1i * l1 * K * 3.14) .* x1,x2),x3);

max_error = max(abs(auto_sol(3:end-1) - exact_sol));
if max_error > 10^-10
    error('Test 4 - 3D to 1D')
end


string_vector_field = 'dot x1 + l1 * Delay(x1,3.14)*x2^2';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

auto_sol = apply(vector_field, x);
auto_sol = verifyfft_in(auto_sol, 1).';
% STILL RETURNS THE IFFT - but I don't seem to be able to update this... whyyyyy??

K = (-2:2);
big_K = -6:6;
big_x1 = [0*(-6:-3), x1, 0*(3:6)];
exact_sol =1i*big_K.*big_x1 + l1 * conv(conv(exp(1i * l1 * K * 3.14) .* x1,x2),x2);

max_error = max(abs(auto_sol(3:end-1) - exact_sol));
if max_error > 10^-10
    error('Test 5 - 3D to 1D')
end


string_vector_field = 'dot x1 + l1 * Delay(x1,3.14)*x1*x2 + 2 * l1 * x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

auto_sol = apply(vector_field, x);
auto_sol = verifyfft_in(auto_sol, 1).';
% STILL RETURNS THE IFFT - but I don't seem to be able to update this... whyyyyy??

K = (-2:2);
exact_sol =1i*K.*x1 + l1 * conv(conv(exp(1i * l1 * K * 3.14) .* x1,x1,'same'),x2,'same') + 2 * l1 * x1;

max_error = max(abs(auto_sol(7:11) - exact_sol));
if max_error > 10^-10
    error('Test 6 - 3D to 1D')
end

disp('The evaluation of a 1D vector field of a 3D vector with delay and convolutions is successful')


%% PART 2: higher dimensions vector field
%
% subpart 2: higher dim vector, higher dim vector field


x1 = [0,-1i,1,+1i,0];
x2 = [1,2,3,2,1] + [0,1,0,-1,0]*1i;
x3 = [2,-1,0,-1,2] + [1,1,0,-1,-1]*1i;

l1 = 1;
x = Xi_vector(l1, [x1;x2;x3]);
string_vector_field = 'dot x1 + l1 * x1 + x2 + x3 \n dot x2 + x3 - x1 \n dot x3';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

auto_sol = apply(vector_field, x);
for i=1:3
    auto_sol(i,:) = verifyfft_in(auto_sol(i,:), 1).';
end
% STILL RETURNS THE IFFT - but I don't seem to be able to update this... whyyyyy??

K = (-2:2);
exact_sol_1 = 1i*K.*x1 + l1 * x1 + x2 + x3;
exact_sol_2 = 1i*K.*x2 + x3 - x1;
exact_sol_3 = 1i*K.*x3;

max_error = max(abs(auto_sol(1,3:end-1) - exact_sol_1)) + ...
    max(abs(auto_sol(2,3:end-1) - exact_sol_2))+ ...
    max(abs(auto_sol(3,3:end-1) - exact_sol_3));
if max_error > 10^-10
    error('Test 1 - 3D to 3D')
end


string_vector_field = 'dot x1 + l1 * x1 * x2 + x2 * x3 \n dot x2 + x3 * x1 - x1 \n dot x3 + x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

auto_sol = apply(vector_field, x);
for i=1:3
    auto_sol(i,:) = verifyfft_in(auto_sol(i,:), 1).';
end
% STILL RETURNS THE IFFT - but I don't seem to be able to update this... whyyyyy??

K = (-2:2);
exact_sol_1 = 1i*K.*x1 + l1 * conv(x1,x2,'same') + conv(x2,x3,'same');
exact_sol_2 = 1i*K.*x2 + conv(x1,x3,'same') - x1;
exact_sol_3 = 1i*K.*x3 + x1 ;

max_error = max(abs(auto_sol(1,7:11) - exact_sol_1)) + ...
    max(abs(auto_sol(2,7:11) - exact_sol_2))+ ...
    max(abs(auto_sol(3,7:11) - exact_sol_3));
if max_error > 10^-10
    error('Test 2 - 3D to 3D')
end


string_vector_field = 'dot x1 + l1 * x1 * x2 + x2^2 * x3 \n dot x2 + x3 * x1^2 - x1 \n dot x3 + x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

auto_sol = apply(vector_field, x);
for i=1:3
    auto_sol(i,:) = verifyfft_in(auto_sol(i,:), 1).';
end
% STILL RETURNS THE IFFT - but I don't seem to be able to update this... whyyyyy??

K = (-2:2);
big_K = -6:6;
big_x1 = [0*(-6:-3), x1, 0*(3:6)];
big_x2 = [0*(-6:-3), x2, 0*(3:6)];
one_vec = [0,0,1,0,0];
exact_sol_1 = 1i*big_K.*big_x1 + l1 * conv(conv(x1,x2),one_vec) + conv(x2,conv(x2,x3));
exact_sol_2 = 1i*big_K.*big_x2 + conv(x1,conv(x1,x3)) - big_x1;
exact_sol_3 = 1i*K.*x3 + x1 ;

max_error = max(abs(auto_sol(1,3:end-1) - exact_sol_1)) + ...
    max(abs(auto_sol(2,3:end-1) - exact_sol_2))+ ...
    max(abs(auto_sol(3,7:11) - exact_sol_3));
if max_error > 10^-10
    error('Test 2b - 3D to 3D')
end


string_vector_field = 'dot x1 + l1 * Delay(x1,3.14) * x1 + x2^2 * x3 \n dot x2 + Delay(x2,3) * x1^2 - x1 \n dot x3 + x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

auto_sol = apply(vector_field, x);
for i=1:3
    auto_sol(i,:) = verifyfft_in(auto_sol(i,:), 1).';
end
% STILL RETURNS THE IFFT - but I don't seem to be able to update this... whyyyyy??

K = (-2:2);
big_K = -6:6;
big_x1 = [0*(-6:-3), x1, 0*(3:6)];
big_x2 = [0*(-6:-3), x2, 0*(3:6)];
one_vec = [0,0,1,0,0];

exact_sol_1 = 1i*big_K.*big_x1 + l1 * conv(conv(exp(1i * l1 * K * 3.14) .* x1,x1),one_vec) + conv(x2,conv(x2,x3));
exact_sol_2 = 1i*big_K.*big_x2 + conv(x1,conv(x1,exp(1i * l1 * K * 3) .* x2)) - big_x1;
exact_sol_3 = 1i*K.*x3 + x1 ;

max_error = max(abs(auto_sol(1,3:end-1) - exact_sol_1)) + ...
    max(abs(auto_sol(2,3:end-1) - exact_sol_2))+ ...
    max(abs(auto_sol(3,7:11) - exact_sol_3));
if max_error > 10^-10
    error('Test 3 - 3D to 3D')
end



string_vector_field = 'dot x1 + l1 * Delay(x1,3.14) * Delay(x2,4) + x2^2 * x3 \n dot x2 + Delay(x2,3) * Delay(x1,2) *x3 - x1 \n dot x3 + x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

auto_sol = apply(vector_field, x);
for i=1:3
    auto_sol(i,:) = verifyfft_in(auto_sol(i,:), 1).';
end

K = (-2:2);
big_K = -6:6;
big_x1 = [0*(-6:-3), x1, 0*(3:6)];
big_x2 = [0*(-6:-3), x2, 0*(3:6)];
one_vec = [0,0,1,0,0];

exact_sol_1 = 1i*big_K.*big_x1 + l1 * conv(conv(exp(1i * l1 * K * 3.14) .* x1,exp(1i * l1 * K * 4) .*x2),one_vec) + conv(x2,conv(x2,x3));
exact_sol_2 = 1i*big_K.*big_x2 + conv(x3,conv(exp(1i * l1 * K * 2) .*x1,exp(1i * l1 * K * 3) .* x2)) - big_x1;
exact_sol_3 = 1i*K.*x3 + x1 ;

max_error = max(abs(auto_sol(1,3:end-1) - exact_sol_1)) + ...
    max(abs(auto_sol(2,3:end-1) - exact_sol_2))+ ...
    max(abs(auto_sol(3,7:11) - exact_sol_3));
if max_error > 10^-10
    error('Test 4 - 3D to 1D')
end

disp('The evaluation of a 3D vector field of a 3D vector with delay and convolutions is successful')


%% PART 3:  derivatives 1D vector field in 1 variable


x1 = [0,-1i,1,+1i,0];
l1 = 1;
x = Xi_vector(l1, x1);
scalar_equation = fancy_scalar_condition(x);

string_vector_field = 'dot x1 + l1 * x1^2 ';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,1);
full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda).' - conv(x1,x1));
max_error(2) = norm(auto_der.derivative_Fx_diagonal - 1);
max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix).' - 2*l1*conv(x1,one_vec));
if ~isempty(auto_der.derivative_Fx_delay{1}{1})
    error('wrong size delay derivative')
elseif ~isempty(auto_der.derivative_Fx_delay{1}{2})
    error('wrong size delay derivative')
end

max_error = norm(max_error);
if max_error > 10^-10
    error('Test 1 delay - 1D')
end


x1 = [0,-1i,1,+1i,0];
l1 = 2;
x = Xi_vector(l1, x1);
string_vector_field = 'l1 * dot x1 + x1^2 + 2 * x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,1);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda).' - conv(1i*K.*x1,one_vec));
max_error(2) = norm(auto_der.derivative_Fx_diagonal - l1);
max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix).' - 2*conv(x1,one_vec) - 2*conv(one_vec,one_vec));
if ~isempty(auto_der.derivative_Fx_delay{1}{1})
    error('wrong size delay derivative')
elseif ~isempty(auto_der.derivative_Fx_delay{1}{2})
    error('wrong size delay derivative')
end

max_error = norm(max_error);
if max_error > 10^-10
    error('Test 2 delay - 1D')
end



x1 = [0,-1i,1,+1i,0];
l1 = 1.5;
x = Xi_vector(l1, x1);
Delay = @(x, tau) exp(1i*l1*(-(length(x)-1)/2:(length(x)-1)/2)*tau).*x;

string_vector_field = 'dot x1 + l1 * Delay(x1,3.14)*x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,1);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda).' - conv(exp(1i*K*l1*3.14).*x1,x1)...
    - l1 * conv(1i*K*3.14.*exp(1i*K*l1*3.14).*x1,x1));
max_error(2) = norm(auto_der.derivative_Fx_diagonal - 1);
max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix).' - l1 * conv(Delay(x1,3.14),one_vec));
if ~isempty(auto_der.derivative_Fx_delay{1}{1})
    error('wrong size delay derivative')
elseif isempty(auto_der.derivative_Fx_delay{1}{2})
    error('wrong size delay derivative')
end
max_error(4) = norm(auto_der.derivative_Fx_delay{1}{2}.convolution - l1*conv(x1,one_vec).');

max_error = norm(max_error);
if max_error > 10^-10
    error('Test 3 delay - 1D')
end



x1 = [0,-1i,1,+1i,0];
l1 = 2;
x = Xi_vector(l1, x1);
Delay = @(x, tau) exp(1i*l1*(-(length(x)-1)/2:(length(x)-1)/2)*tau).*x;

string_vector_field = 'l1 * Delay(x1,3.14)*x1^2';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,1);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda).' - conv(Delay(x1,3.14),conv(x1,x1))...
    - l1 * conv((1i*K*3.14).*Delay(x1,3.14),conv(x1,x1)) );
max_error(2) = norm(auto_der.derivative_Fx_diagonal - 0);
max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix).' - l1 * 2 * conv(conv(Delay(x1,3.14),x1),one_vec));
if isempty(auto_der.derivative_Fx_delay{1}{1})
    error('wrong size delay derivative')
end
max_error(4) = norm(auto_der.derivative_Fx_delay{1}{1}.convolution - l1*conv(conv(x1,x1),one_vec).');
max_error(5) = norm(auto_der.derivative_Fx_delay{1}{1}.exp_coef - 1i*3.14*l1);

if max(max_error) > 10^-10
    error('Test 4 delay derivative - 1D')
end

x1 = [0,-1i,1,+1i,0];
l1 = 1.5;
x = Xi_vector(l1, x1);
Delay = @(x, tau) exp(1i*l1*(-(length(x)-1)/2:(length(x)-1)/2)*tau).*x;
string_vector_field = 'dot x1 + l1 * Delay(x1,3.14)*x1^2';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,1);
full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda).' - conv(Delay(x1,3.14),conv(x1,x1))...
    - l1 * conv((1i*K*3.14).*Delay(x1,3.14),conv(x1,x1)));
max_error(2) = norm(auto_der.derivative_Fx_diagonal - 1);
max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix).' - l1 * 2 * conv(conv(Delay(x1,3.14),x1),one_vec));
if ~isempty(auto_der.derivative_Fx_delay{1}{1})
    error('wrong size delay derivative')
end
if isempty(auto_der.derivative_Fx_delay{1}{2})
    error('wrong size delay derivative')
end
max_error(4) = norm(auto_der.derivative_Fx_delay{1}{2}.convolution - l1*conv(conv(x1,x1),one_vec).');
max_error(5) = norm(auto_der.derivative_Fx_delay{1}{2}.exp_coef - 1i*3.14*l1);


if max(max_error) > 10^-10
    error('Test 5 derivative delay - 1D')
end

x1 = [0,-1i,1,+1i,0];
l1 = 1.5;
x = Xi_vector(l1, x1);
Delay = @(x, tau) exp(1i*l1*(-(length(x)-1)/2:(length(x)-1)/2)*tau).*x;

string_vector_field = 'l1* dot x1 + l1 * Delay(x1,3.34)*x1^2 + 2 * l1 * x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,1);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda).' - ...
    1i*conv(K.*x1, conv(one_vec,one_vec)) - conv(Delay(x1,3.34),conv(x1,x1)) -  2 * conv(x1, conv(one_vec,one_vec))...
    - l1 * conv((1i*K*3.34).*Delay(x1,3.34),conv(x1,x1)));
max_error(2) = norm(auto_der.derivative_Fx_diagonal - l1);
max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix).' - ...
    l1 * 2 * conv(conv(Delay(x1,3.34),x1),one_vec) - 2 * l1 * conv(one_vec, conv(one_vec,one_vec)));

if ~isempty(auto_der.derivative_Fx_delay{1}{1})
    error('wrong size delay derivative')
end
if ~isempty(auto_der.derivative_Fx_delay{1}{3})
    error('wrong size delay derivative')
end
if isempty(auto_der.derivative_Fx_delay{1}{2})
    error('wrong size delay derivative')
end
max_error(4) = norm(auto_der.derivative_Fx_delay{1}{2}.convolution - l1*conv(x1, conv(x1,one_vec)).');
max_error(5) = norm(auto_der.derivative_Fx_delay{1}{2}.exp_coef - 1i*3.34*l1);

if max(max_error) > 10^-10
    error('Test 6 delay derivative - 1D')
end

disp('The evaluation of the derivative of a 1D vector field with delay and convolutions is successful')


%% PART 4: derivatives of higher dimensions vector field
%
% subpart 1: higher dim vector, 1D vector field


x1 = [0,-1i,1,+1i,0];
x2 = [1,2,3,2,1] + [0,1,0,-1,0]*1i;
x3 = [2,-1,0,-1,2] + [1,1,0,-1,-1]*1i;

l1 = -2;
x = Xi_vector(l1, [x1;x2;x3]);
Delay = @(x, tau) exp(1i*l1*(-(length(x)-1)/2:(length(x)-1)/2)*tau).*x;

scalar_equation = fancy_scalar_condition(x);

string_vector_field = 'dot x1 + l1 * x1 + x2 + x3 ';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda).' - x1);
max_error(2) = norm(auto_der.derivative_Fx_diagonal - [1;0;0]);
max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix(1,:,:)).' - l1*one_vec) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,:,:)).' - one_vec) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,:,:)).' - one_vec);

if ~isempty(auto_der.derivative_Fx_delay{1}{1})
    error('wrong size delay derivative')
end
if ~isempty(auto_der.derivative_Fx_delay{1}{2})
    error('wrong size delay derivative')
end
if ~isempty(auto_der.derivative_Fx_delay{1}{3})
    error('wrong size delay derivative')
end

if max(max_error) > 10^-10
    error('Test 1 delay derivative - 3D to 1D')
end


string_vector_field = 'dot x2 + l1 * x1^2 + 2 * x1 * x2';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda).' - conv(x1,x1));
max_error(2) = norm(auto_der.derivative_Fx_diagonal - [0;1;0]);
max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix(1,:,:)).' - l1*2*conv(one_vec,x1) - 2*conv(one_vec,x2) ) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,:,:)).' - 2*conv(x1,one_vec)) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,:,:)).' - 0);

if ~isempty(auto_der.derivative_Fx_delay{1}{1})
    error('wrong size delay derivative')
end
if ~isempty(auto_der.derivative_Fx_delay{1}{2})
    error('wrong size delay derivative')
end
if ~isempty(auto_der.derivative_Fx_delay{1}{3})
    error('wrong size delay derivative')
end

if max(max_error) > 10^-10
    error('Test 2 delay derivative - 3D to 1D')
end


string_vector_field = 'dot x1 + l1 * Delay(x1,3.14)*x2';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda).' - conv(Delay(x1,3.14),x2)...
    - l1 * conv((1i*K*3.14).*Delay(x1,3.14),x2));
max_error(2) = norm(auto_der.derivative_Fx_diagonal - [1;0;0]);
max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix(1,:,:)).' - 0) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,:,:)).' - l1*conv(Delay(x1,3.14),one_vec)) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,:,:)).' - 0);

if ~isempty(auto_der.derivative_Fx_delay{1}{1})
    error('wrong size delay derivative')
end
if isempty(auto_der.derivative_Fx_delay{1}{2})
    error('wrong size delay derivative')
end
max_error(4) = norm(auto_der.derivative_Fx_delay{1}{2}.convolution - l1*conv(x2,one_vec).');
max_error(5) = norm(auto_der.derivative_Fx_delay{1}{2}.exp_coef - 1i*3.14*l1);

if max(max_error) > 10^-10
    error('Test 3 derivative delay - 3D to 1D')
end



string_vector_field = 'l1 * Delay(x1,3.14)*x2*x3';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda).' - conv(conv(Delay(x1,3.14),x2),x3)...
    - l1 * conv((1i*K*3.14).*Delay(x1,3.14),conv(x2,x3)));
max_error(2) = norm(auto_der.derivative_Fx_diagonal - [0;0;0]);
max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix(1,:,:)).' - 0) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,:,:)).' - l1*conv(conv(Delay(x1,3.14),x3),one_vec)) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,:,:)).' - l1*conv(conv(Delay(x1,3.14),x2),one_vec));

if isempty(auto_der.derivative_Fx_delay{1}{1})
    error('wrong size delay derivative')
end
max_error(4) = norm(auto_der.derivative_Fx_delay{1}{1}.convolution - l1*conv(conv(x2,x3),one_vec).');
max_error(5) = norm(auto_der.derivative_Fx_delay{1}{1}.exp_coef - 1i*3.14*l1);

if max(max_error) > 10^-10
    error('Test 4 - 3D to 1D')
end


string_vector_field = 'l1 * dot x1 + l1 * Delay(x1,3.74)*x2^2';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda).' ...
    - 1i*conv(conv(K.*x1, one_vec),one_vec) - conv(conv(Delay(x1,3.74),x2),x2)...
    - l1 * conv((1i*K*3.74).*Delay(x1,3.74),conv(x2,x2)));
max_error(2) = norm(auto_der.derivative_Fx_diagonal - [l1;0;0]);
max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix(1,:,:)).' - 0) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,:,:)).' - 2*l1*conv(conv(Delay(x1,3.74),x2),one_vec)) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,:,:)).' - 0);

if ~isempty(auto_der.derivative_Fx_delay{1}{1})
    error('wrong size delay derivative')
end
if isempty(auto_der.derivative_Fx_delay{1}{2})
    error('wrong size delay derivative')
end

max_error(4) = norm(auto_der.derivative_Fx_delay{1}{2}.convolution - l1*conv(conv(x2,x2),one_vec).');
max_error(5) = norm(auto_der.derivative_Fx_delay{1}{2}.exp_coef - 1i*3.74*l1);

if max(max_error) > 10^-10
    error('Test 5 delay derivative - 3D to 1D')
end


string_vector_field = 'dot x1 + l1 * Delay(x1,3.14)*x1*x2 + 2 * l1 * x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);
full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda).' ...
    - conv(conv(Delay(x1,3.14), x1),x2) - 2*conv(conv(x1,one_vec),one_vec)...
    - l1 * conv((1i*K*3.14).*Delay(x1,3.14),conv(x1,x2)));
max_error(2) = norm(auto_der.derivative_Fx_diagonal - [1;0;0]);
max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix(1,:,:)).' - l1*conv(conv(Delay(x1,3.14), one_vec),x2) - 2*l1*conv(conv(one_vec,one_vec),one_vec)) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,:,:)).' - l1*conv(conv(Delay(x1,3.14), one_vec),x1)) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,:,:)).' - 0);

if ~isempty(auto_der.derivative_Fx_delay{1}{1})
    error('wrong size delay derivative')
end
if ~isempty(auto_der.derivative_Fx_delay{1}{3})
    error('wrong size delay derivative')
end
if isempty(auto_der.derivative_Fx_delay{1}{2})
    error('wrong size delay derivative')
end

max_error(4) = norm(auto_der.derivative_Fx_delay{1}{2}.convolution - l1*conv(conv(x1,x2),one_vec).');
max_error(5) = norm(auto_der.derivative_Fx_delay{1}{2}.exp_coef - 1i*3.14*l1);

if max(max_error) > 10^-10
    error('Test 6 derivative delay - 3D to 1D')
end

disp('The evaluation of a 1D vector field of a 3D vector with delay and convolutions is successful')


%% PART 4: derivatives of higher dimensions vector field
%
% subpart 2: higher dim vector, higher dim vector field


x1 = [0,-1i,1,+1i,0];
x2 = [1,2,3,2,1] + [0,1,0,-1,0]*1i;
x3 = [2,-1,0,-1,2] + [1,1,0,-1,-1]*1i;

l1 = 1.8;
x = Xi_vector(l1, [x1;x2;x3]);
Delay = @(x, tau) exp(1i*l1*(-(length(x)-1)/2:(length(x)-1)/2)*tau).*x;
string_vector_field = 'dot x1 + l1 * x1 + x2 + x3 \n l1 * dot x2 + x3 - x1 \n dot x3';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);
full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda(:,1,:)).' - x1)+...
    norm(squeeze(auto_der.derivative_Flambda(:,2,:)).' - 1i*K.*x2)+...
    norm(squeeze(auto_der.derivative_Flambda(:,3,:)).' - 0);

max_error(2) = norm(auto_der.derivative_Fx_diagonal - [1,0,0;0,l1,0;0,0,1]);
max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix(1,1,:)).' - l1* one_vec) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,1,:)).' - one_vec) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,1,:)).' - one_vec);

max_error(3) = max_error(3) + norm(squeeze(auto_der.derivative_Fx_toeplix(1,2,:)).' - (-1*one_vec)) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,2,:)).' - 0) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,2,:)).' - one_vec);

max_error(3) = max_error(3) + norm(squeeze(auto_der.derivative_Fx_toeplix(1,3,:)).' - 0) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,3,:)).' - 0) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,3,:)).' - 0);

if ~isempty(auto_der.derivative_Fx_delay{1}{1})
    error('wrong size delay derivative')
end
if ~isempty(auto_der.derivative_Fx_delay{1}{3})
    error('wrong size delay derivative')
end
if ~isempty(auto_der.derivative_Fx_delay{1}{2})
    error('wrong size delay derivative')
end

max_error(4) = 0;
max_error(5) = 0;

if max(max_error) > 10^-10
    error('Test 1 derivative delay - 3D to 3D')
end


string_vector_field = 'dot x1 + l1 * x1 * x2 + x2 * x3 \n -4 *dot x2 + x3 * x1 - x1 \n dot x3 + x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda(:,1,:)).' - conv(x1,x2))+...
    norm(squeeze(auto_der.derivative_Flambda(:,2,:)).' - 0)+...
    norm(squeeze(auto_der.derivative_Flambda(:,3,:)).' - 0);

max_error(2) = norm(auto_der.derivative_Fx_diagonal - [1,0,0;0,-4,0;0,0,1]);
max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix(1,1,:)).' - l1* conv(x2,one_vec)) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,1,:)).' - l1* conv(x1,one_vec) - conv(x3,one_vec)) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,1,:)).' - conv(x2,one_vec));

max_error(3) = max_error(3) + norm(squeeze(auto_der.derivative_Fx_toeplix(1,2,:)).' - conv(x3,one_vec) + conv(one_vec,one_vec)) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,2,:)).' - 0) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,2,:)).' - conv(x1,one_vec));

max_error(3) = max_error(3) + norm(squeeze(auto_der.derivative_Fx_toeplix(1,3,:)).' - conv(one_vec,one_vec)) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,3,:)).' - 0) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,3,:)).' - 0);


for i = 1:3
    for j = 1:3
        for k = 1:3
            if i == 3 && k == 3
                continue
            end
            if ~isempty(auto_der.derivative_Fx_delay{i,j}{k})
                error('wrong size delay derivative')
            end
        end
    end
end

max_error(4) = 0;
max_error(5) = 0;

if max_error > 10^-10
    error('Test 2 delay derivative - 3D to 3D')
end



string_vector_field = 'dot x1 + l1 * x1 * x2 + x2^2 * x3 \n -4*dot x2 + x3 * x1^2 - x1 \n dot x3 + x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda(:,1,:)).' - conv(conv(x1,x2),one_vec))+...
    norm(squeeze(auto_der.derivative_Flambda(:,2,:)).' - 0)+...
    norm(squeeze(auto_der.derivative_Flambda(:,3,:)).' - 0);

max_error(2) = norm(auto_der.derivative_Fx_diagonal - [1,0,0;0,-4,0;0,0,1]);

max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix(1,1,:)).' - l1* conv(conv(x2,one_vec),one_vec)) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,1,:)).' - l1* conv(conv(x1,one_vec),one_vec) - 2*conv(conv(x3,x2),one_vec)) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,1,:)).' - conv(conv(x2,x2),one_vec));

max_error(3) = max_error(3) + norm(squeeze(auto_der.derivative_Fx_toeplix(1,2,:)).' - 2*conv(conv(x3,one_vec),x1) + conv(one_vec,conv(one_vec,one_vec))) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,2,:)).' - 0) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,2,:)).' - conv(x1,conv(x1,one_vec)));

max_error(3) = max_error(3) + norm(squeeze(auto_der.derivative_Fx_toeplix(1,3,:)).' - conv(one_vec,conv(one_vec,one_vec))) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,3,:)).' - 0) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,3,:)).' - 0);


for i = 1:3
    for j = 1:3
        for k = 1:3
            if i == 3 && k == 3
                continue
            end
            if ~isempty(auto_der.derivative_Fx_delay{i,j}{k})
                error('wrong size delay derivative')
            end
        end
    end
end

max_error(4) = 0;
max_error(5) = 0;

if max(max_error) > 10^-10
    error('Test 2b derivative delay - 3D to 3D')
end


string_vector_field = 'dot x1 + l1 * Delay(x1,3.14) * x1 + x2^2 * x3 \n dot x2 + Delay(x2,3) * x1^2 - x1 \n l1 * dot x3 + x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda(:,1,:)).' - conv(conv(Delay(x1,3.14), x1),one_vec)- l1*conv((1i*K*3.14).*Delay(x1,3.14),conv(x1,one_vec)))+...
    norm(squeeze(auto_der.derivative_Flambda(:,2,:)).' - conv((1i*K*3).*Delay(x2,3),conv(x1,x1)))+...
    norm(squeeze(auto_der.derivative_Flambda(:,3,:)).' - conv(conv(1i*K.*x3, one_vec),one_vec)...
    );

max_error(2) = norm(auto_der.derivative_Fx_diagonal - [1,0,0;0,1,0;0,0,l1]);

max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix(1,1,:)).' - l1* conv(conv(Delay(x1,3.14),one_vec),one_vec)) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,1,:)).' - 2 * conv(conv(x2,x3),one_vec)) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,1,:)).' - conv(conv(x2,x2),one_vec));

max_error(3) = max_error(3) + norm(squeeze(auto_der.derivative_Fx_toeplix(1,2,:)).' - 2*conv(conv(Delay(x2,3),x1),one_vec) + conv(one_vec, conv(one_vec, one_vec))) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,2,:)).' - 0) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,2,:)).' - 0);

max_error(3) = max_error(3) + norm(squeeze(auto_der.derivative_Fx_toeplix(1,3,:)).' - conv(one_vec,conv(one_vec,one_vec))) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,3,:)).' - 0) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,3,:)).' - 0);


for i = 1:3
    for j = 1:3
        for k = 1:3
            if i == 3 && k == 3
                continue
            end
            if all([i,j,k] == [1,1,2]) || all([i,j,k] == [2,2,2])
                if isempty(auto_der.derivative_Fx_delay{i,j}{k})
                    disp([i,j,k])
                    error('wrong size delay derivative')
                end
            else
                if ~isempty(auto_der.derivative_Fx_delay{i,j}{k})
                    disp([i,j,k])
                    error('wrong size delay derivative')
                end
            end
        end
    end
end
% term 1,1,2 : l1 * Delay(x1,3.14) * x1
% term 2,2,2 : Delay(x2,3) * x1^2
max_error(4) = norm(auto_der.derivative_Fx_delay{1,1}{2}.convolution.' - l1* conv(x1, conv(one_vec,one_vec))) + ...
    norm(auto_der.derivative_Fx_delay{2,2}{2}.convolution.' - conv(x1, conv(x1,one_vec)));

max_error(5) = norm(auto_der.derivative_Fx_delay{1,1}{2}.exp_coef - 1i*3.14*l1) + ...
    norm(auto_der.derivative_Fx_delay{2,2}{2}.exp_coef - 1i*3*l1);

if max(max_error) > 10^-10
    error('Test 3 delay derivative - 3D to 3D')
end

% dot x1 + l1 * Delay(x1,3.14) * Delay(x2,4) + x2^2 * x3
% dot x2 + Delay(x2,3) * Delay(x1,2) * x3 - x1 
% dot x3 + x1
string_vector_field = 'dot x1 + l1 * Delay(x1,3.14) * Delay(x2,4) + x2^2 * x3 \n dot x2 + Delay(x2,3) * Delay(x1,2) *x3 - x1 \n dot x3 + x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda(1,1,:)).' ...
    - conv(conv(Delay(x1,3.14), Delay(x2,4)),one_vec) ...
    - l1 * conv((1i*K*3.14).*Delay(x1,3.14),conv(Delay(x2,4),one_vec))...
    - l1 * conv(Delay(x1,3.14),conv((1i*K*4).*Delay(x2,4),one_vec)))+...
    norm(squeeze(auto_der.derivative_Flambda(1,2,:)).'...
    - conv((1i*K*2).*Delay(x1,2),conv(Delay(x2,3),x3))...
    - conv(Delay(x1,2),conv((1i*K*3).*Delay(x2,3),x3)) )+...
    norm(squeeze(auto_der.derivative_Flambda(1,3,:)).' - 0);

max_error(2) = norm(auto_der.derivative_Fx_diagonal - [1,0,0;0,1,0;0,0,1]);

max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix(1,1,:)).' - 0) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,1,:)).' - 2 * conv(conv(x2,x3),one_vec)) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,1,:)).' - conv(conv(x2,x2),one_vec));

max_error(3) = max_error(3) + norm(squeeze(auto_der.derivative_Fx_toeplix(1,2,:)).' + conv(one_vec, conv(one_vec, one_vec))) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,2,:)).' - 0) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,2,:)).' - conv(conv(Delay(x2,3), Delay(x1,2)),one_vec));

max_error(3) = max_error(3) + norm(squeeze(auto_der.derivative_Fx_toeplix(1,3,:)).' - conv(one_vec,conv(one_vec,one_vec))) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,3,:)).' - 0) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,3,:)).' - 0);


for i = 1:3
    for j = 1:3
        for k = 1:3
            if i == 3 && k == 3
                continue
            end
            if all([i,j,k] == [1,1,2]) || all([i,j,k] == [1,2,2]) || all([i,j,k] == [2,1,2]) || all([i,j,k] == [2,2,2])
                if isempty(auto_der.derivative_Fx_delay{i,j}{k})
                    disp([i,j,k])
                    error('wrong size delay derivative')
                end
            else
                if ~isempty(auto_der.derivative_Fx_delay{i,j}{k})
                    disp([i,j,k])
                    error('wrong size delay derivative')
                end
            end
        end
    end
end
% term 1,1,2 : l1 * Delay(x1,3.14) * Delay(x2,4)
% term 2,2,2 : Delay(x2,3) * Delay(x1,2) *x3
max_error(4) = norm(auto_der.derivative_Fx_delay{1,1}{2}.convolution.' - l1* conv(Delay(x2,4), conv(one_vec,one_vec))) + ...
    norm(auto_der.derivative_Fx_delay{1,2}{2}.convolution.' - l1 * conv(Delay(x1,3.14), conv(one_vec,one_vec))) + ...
    norm(auto_der.derivative_Fx_delay{2,1}{2}.convolution.' - conv(Delay(x2,3), conv(x3,one_vec))) + ...
    norm(auto_der.derivative_Fx_delay{2,2}{2}.convolution.' - conv(Delay(x1,2), conv(x3,one_vec)));


max_error(5) = norm(auto_der.derivative_Fx_delay{1,1}{2}.exp_coef - 1i*3.14*l1) + ...
    norm(auto_der.derivative_Fx_delay{1,2}{2}.exp_coef - 1i*4*l1) + ...
    norm(auto_der.derivative_Fx_delay{2,1}{2}.exp_coef - 1i*2*l1) + ...
    norm(auto_der.derivative_Fx_delay{2,2}{2}.exp_coef - 1i*3*l1);
if max(max_error) > 10^-10
    error('Test 4 derivative delay - 3D to 3D')
end

disp('The evaluation of a 3D vector field of a 3D vector with delay and convolutions is successful')


%% PART 5: composing derivatives 1D vector field in 1 variable

conv_mat = @(x) conj(toeplitz(x((length(x)+1)/2:end)));

conv_mat_big =@(x) conj(toeplitz([x((length(x)+1)/2:end),0*x((length(x)+3)/2:end)]));


x1 = [0,-1i,1,+1i,0];
l1 = 1;
x = Xi_vector(l1, x1);
scalar_equation = fancy_scalar_condition(x);

string_vector_field = 'dot x1 + l1 * x1^2 ';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,1);
full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);
auto_mat_der = derivative_to_matrix(auto_der);

mat_K = diag(-2:2);
mat_bigK = diag(-4:4);
mat_verybigK = diag(-6:6);

scal_coefs = squeeze(scalar_equation.linear_coef{2}).';
hand_solution = [0, conv(scal_coefs,one_vec);
    conv(x1,x1).', 1i*mat_bigK + 2*l1 * conv_mat_big(conv(x1,one_vec))];

max_error = norm(hand_solution - auto_mat_der);

if max_error > 10^-10
    error('Test 1 delay - 1D')
end


x1 = [0,-1i,1,+1i,0];
l1 = 2;
x = Xi_vector(l1, x1);
string_vector_field = 'l1 * dot x1 + x1^2 + 2 * x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,1);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

auto_mat_der = derivative_to_matrix(auto_der);

scal_coefs = squeeze(scalar_equation.linear_coef{2}).';
hand_solution = [0, conv(scal_coefs,one_vec);
    conv(1i*K.*x1,one_vec).', l1*1i*mat_bigK + conv_mat_big(2*conv(x1,one_vec) + 2*conv(one_vec,one_vec))];

max_error = norm(hand_solution - auto_mat_der);

if max_error > 10^-10
    error('Test 2 delay - 1D')
end



x1 = [0,-1i,1,+1i,0];
l1 = 1.5;
x = Xi_vector(l1, x1);
Delay = @(x, tau) exp(1i*l1*(-(length(x)-1)/2:(length(x)-1)/2)*tau).*x;

string_vector_field = 'Delay(x1,3.14)';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,1);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda).' - 1i * K * 3.14 .* Delay(x1,3.14));
max_error(2) = norm(auto_der.derivative_Fx_diagonal - 0);
max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix).' - 0);
max_error(4) = norm(auto_der.derivative_Fx_delay{1}{1}.convolution - one_vec.');
max_error(5) = norm(auto_der.derivative_Fx_delay{1}{1}.exp_coef - 1i*3.14*l1);

auto_mat_der = derivative_to_matrix(auto_der);

scal_coefs = squeeze(scalar_equation.linear_coef{2}).';
delay_mat = expm(1i*3.14*l1 * mat_K);
hand_solution = [0, scal_coefs;
    (1i * K * 3.14 .* Delay(x1,3.14)).', delay_mat];

max_error = norm(hand_solution - auto_mat_der);

max_error = norm(max_error);
if max_error > 10^-10
    error('Test 3 delay - 1D')
end


x1 = [0,-1i,1,+1i,0];
l1 = 2;
x = Xi_vector(l1, x1);
Delay = @(x, tau) exp(1i*l1*(-(length(x)-1)/2:(length(x)-1)/2)*tau).*x;

string_vector_field = 'l1 * Delay(x1,3.14)*x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,1);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda).' - conv(Delay(x1,3.14),x1) - l1 * conv(1i *K * 3.14 .*Delay(x1,3.14),x1));
max_error(2) = norm(auto_der.derivative_Fx_diagonal - 0);
max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix).' - l1 * conv(Delay(x1,3.14),one_vec));
if isempty(auto_der.derivative_Fx_delay{1}{1})
    error('wrong size delay derivative')
end
max_error(4) = norm(auto_der.derivative_Fx_delay{1}{1}.convolution - l1*conv(one_vec,x1).');
max_error(5) = norm(auto_der.derivative_Fx_delay{1}{1}.exp_coef - 1i*3.14*l1);

auto_mat_der = derivative_to_matrix(auto_der);

scal_coefs = squeeze(scalar_equation.linear_coef{2}).';
delay_mat = conv_mat_big(l1*conv(x1,one_vec)) *  expm(1i*3.14*l1 * mat_bigK);
dot_der = 0;
mat_conv = conv_mat_big(l1 * conv(Delay(x1,3.14),one_vec));

der_lambda = conv(Delay(x1,3.14),x1) + l1 * conv(1i *K * 3.14 .*Delay(x1,3.14),x1);
hand_solution = [0, conv(scal_coefs,one_vec);
    der_lambda.', dot_der + delay_mat+mat_conv];

max_error = norm(hand_solution - auto_mat_der);


if max(max_error) > 10^-10
    error('Test 4 delay derivative - 1D')
end




x1 = [0,-1i,1,+1i,0];
l1 = 2;
x = Xi_vector(l1, x1);
Delay = @(x, tau) exp(1i*l1*(-(length(x)-1)/2:(length(x)-1)/2)*tau).*x;

string_vector_field = '0*l1 * dot x1 - l1 * Delay(x1,3.14)*x1^2';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,1);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda).' - 0*conv(one_vec, conv(one_vec, 1i*K.* x1)) + conv(Delay(x1,3.14),conv(x1,x1)) + l1 * conv(1i *K*3.14.*Delay(x1,3.14),conv(x1,x1)));
max_error(2) = norm(auto_der.derivative_Fx_diagonal - 0);
max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix).' + l1 * 2 * conv(conv(Delay(x1,3.14),x1),one_vec));
if isempty(auto_der.derivative_Fx_delay{1}{2})
    error('wrong size delay derivative')
end
max_error(4) = norm(auto_der.derivative_Fx_delay{1}{2}.convolution + l1*conv(conv(x1,x1),one_vec).');
max_error(5) = norm(auto_der.derivative_Fx_delay{1}{2}.exp_coef - 1i*3.14*l1);

auto_mat_der = derivative_to_matrix(auto_der);

scal_coefs = squeeze(scalar_equation.linear_coef{2}).';
delay_mat = - conv_mat_big(l1*conv(conv(x1,x1),one_vec)) *  expm(1i*3.14*l1 * mat_verybigK);
dot_der = 0*mat_verybigK * 1i * l1;
mat_conv = - conv_mat_big(l1 * 2 * conv(conv(Delay(x1,3.14),x1),one_vec));

d_lambda = - conv(Delay(x1,3.14),conv(x1,x1)) -  l1 * conv(1i *K*3.14.*Delay(x1,3.14),conv(x1,x1));
hand_solution = [0, conv(one_vec,conv(scal_coefs,one_vec));
    d_lambda.', dot_der + delay_mat+mat_conv];

max_error = norm(hand_solution - auto_mat_der);

if max(max_error) > 10^-10
    error('Test 4 delay derivative - 1D')
end

x1 = [0,-1i,1,+1i,0];
l1 = 2;
x = Xi_vector(l1, x1);
Delay = @(x, tau) exp(1i*l1*(-(length(x)-1)/2:(length(x)-1)/2)*tau).*x;

string_vector_field = 'l1 * dot x1 - l1 * Delay(x1,3.14)*x1^2';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,1);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda).' - conv(one_vec, conv(one_vec, 1i*K.* x1)) + conv(Delay(x1,3.14),conv(x1,x1))+ l1 * conv(1i *K*3.14.*Delay(x1,3.14),conv(x1,x1)));
max_error(2) = norm(auto_der.derivative_Fx_diagonal - l1);
max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix).' + l1 * 2 * conv(conv(Delay(x1,3.14),x1),one_vec));
if isempty(auto_der.derivative_Fx_delay{1}{2})
    error('wrong size delay derivative')
end
max_error(4) = norm(auto_der.derivative_Fx_delay{1}{2}.convolution + l1*conv(conv(x1,x1),one_vec).');
max_error(5) = norm(auto_der.derivative_Fx_delay{1}{2}.exp_coef - 1i*3.14*l1);

auto_mat_der = derivative_to_matrix(auto_der);

scal_coefs = squeeze(scalar_equation.linear_coef{2}).';
delay_mat = - conv_mat_big(l1*conv(conv(x1,x1),one_vec)) *  expm(1i*3.14*l1 * mat_verybigK);
dot_der = mat_verybigK * 1i * l1;
mat_conv = - conv_mat_big(l1 * 2 * conv(conv(Delay(x1,3.14),x1),one_vec));

d_lambda = conv(one_vec, conv(one_vec, 1i*K.* x1)) - conv(Delay(x1,3.14),conv(x1,x1))- l1 * conv(1i *K*3.14.*Delay(x1,3.14),conv(x1,x1));
hand_solution = [0, conv(one_vec,conv(scal_coefs,one_vec));
    d_lambda.', dot_der + delay_mat+mat_conv];

max_error = norm(hand_solution - auto_mat_der);

if max(max_error) > 10^-10
    error('Test 5 derivative delay - 1D')
end


x1 = [0,-1i,1,+1i,0];
l1 = 1.5;
x = Xi_vector(l1, x1);
Delay = @(x, tau) exp(1i*l1*(-(length(x)-1)/2:(length(x)-1)/2)*tau).*x;

string_vector_field = 'l1* dot x1 + l1 * Delay(x1,3.34)*x1^2 + 2 * l1 * x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,1);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda).' - ...
    1i*conv(K.*x1, conv(one_vec,one_vec)) - conv(Delay(x1,3.34),conv(x1,x1)) -  2 * conv(x1, conv(one_vec,one_vec))...
    - l1 * conv(1i *K*3.34.*Delay(x1,3.34),conv(x1,x1)));
max_error(2) = norm(auto_der.derivative_Fx_diagonal - l1);
max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix).' - ...
    l1 * 2 * conv(conv(Delay(x1,3.34),x1),one_vec) - 2 * l1 * conv(one_vec, conv(one_vec,one_vec)));

if ~isempty(auto_der.derivative_Fx_delay{1}{1})
    error('wrong size delay derivative')
end
if ~isempty(auto_der.derivative_Fx_delay{1}{3})
    error('wrong size delay derivative')
end
if isempty(auto_der.derivative_Fx_delay{1}{2})
    error('wrong size delay derivative')
end
max_error(4) = norm(auto_der.derivative_Fx_delay{1}{2}.convolution - l1*conv(x1, conv(x1,one_vec)).');
max_error(5) = norm(auto_der.derivative_Fx_delay{1}{2}.exp_coef - 1i*3.34*l1);

auto_mat_der = derivative_to_matrix(auto_der);

scal_coefs = squeeze(scalar_equation.linear_coef{2}).';
delay_mat = conv_mat_big(l1*conv(x1, conv(x1,one_vec))) *  expm(1i*3.34*l1 * mat_verybigK);
dot_der = mat_verybigK * 1i * l1;
mat_conv = conv_mat_big(l1 * 2 * conv(conv(Delay(x1,3.34),x1),one_vec) + 2 * l1 * conv(one_vec, conv(one_vec,one_vec)));

d_lambda = 1i*conv(K.*x1, conv(one_vec,one_vec)) + conv(Delay(x1,3.34),conv(x1,x1)) +  2 * conv(x1, conv(one_vec,one_vec))+l1 * conv(1i *K*3.34.*Delay(x1,3.34),conv(x1,x1));
hand_solution = [0, conv(one_vec,conv(scal_coefs,one_vec));
    d_lambda.', dot_der + delay_mat+mat_conv];

max_error = norm(hand_solution - auto_mat_der);

if max(max_error) > 10^-10
    error('Test 6 delay derivative - 1D')
end

disp('The evaluation of the derivative of a 1D vector field with delay and convolutions is successful')


%% PART 4: derivatives of higher dimensions vector field
%
% subpart 1: higher dim vector, 1D vector field
% can't be tested
%
% PART 4: derivatives of higher dimensions vector field
%
% subpart 2: higher dim vector, higher dim vector field


x1 = [0,-1i,1,+1i,0];
x2 = [1,2,3,2,1] + [0,1,0,-1,0]*1i;
x3 = [2,-1,0,-1,2] + [1,1,0,-1,-1]*1i;

l1 = 1.8;
x = Xi_vector(l1, [x1;x2;x3]);
Delay = @(x, tau) exp(1i*l1*(-(length(x)-1)/2:(length(x)-1)/2)*tau).*x;
string_vector_field = 'dot x1 + l1 * x1 + x2 + x3 \n l1 * dot x2 + x3 - x1 \n dot x3';

scalar_equation = fancy_scalar_condition(x);
vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);
full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda(:,1,:)).' - x1)+...
    norm(squeeze(auto_der.derivative_Flambda(:,2,:)).' - 1i*K.*x2)+...
    norm(squeeze(auto_der.derivative_Flambda(:,3,:)).' - 0);

max_error(2) = norm(auto_der.derivative_Fx_diagonal - [1,0,0;0,l1,0;0,0,1]);
max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix(1,1,:)).' - l1* one_vec) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,1,:)).' - one_vec) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,1,:)).' - one_vec);

max_error(3) = max_error(3) + norm(squeeze(auto_der.derivative_Fx_toeplix(1,2,:)).' - (-1*one_vec)) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,2,:)).' - 0) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,2,:)).' - one_vec);

max_error(3) = max_error(3) + norm(squeeze(auto_der.derivative_Fx_toeplix(1,3,:)).' - 0) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,3,:)).' - 0) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,3,:)).' - 0);

if ~isempty(auto_der.derivative_Fx_delay{1}{1})
    error('wrong size delay derivative')
end
if ~isempty(auto_der.derivative_Fx_delay{1}{3})
    error('wrong size delay derivative')
end
if ~isempty(auto_der.derivative_Fx_delay{1}{2})
    error('wrong size delay derivative')
end

max_error(4) = 0;
max_error(5) = 0;

% string_vector_field = 'dot x1 + l1 * x1 + x2 + x3 \n l1 * dot x2 + x3 - x1 \n dot x3';

auto_mat_der = derivative_to_matrix(auto_der);

scal_coefs1 = squeeze(scalar_equation.linear_coef{2}(:,1,:)).';
scal_coefs2 = squeeze(scalar_equation.linear_coef{2}(:,2,:)).';
scal_coefs3 = squeeze(scalar_equation.linear_coef{2}(:,3,:)).';
scal_coef = [scal_coefs1, scal_coefs2, scal_coefs3];
delay_mat = 0;
dot_der = kron([1,0,0; 0,l1,0;0,0,1], mat_K * 1i);
mat_conv = kron( [l1,1,1;-1,0,1;0,0,0], conv_mat_big(one_vec));
d_lambda = [x1, 1i*K.*x2, 0*x3].';

hand_solution = [0, scal_coef;
    d_lambda, dot_der + mat_conv];

max_error = norm(hand_solution - auto_mat_der);

if max(max_error) > 10^-10
    error('Test 1 derivative delay - 3D to 3D')
end



string_vector_field = 'dot x1 + l1 * x1 * x2 + x2 * x3 \n -4 *dot x2 + x3 * x1 - x1 \n dot x3 + x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda(:,1,:)).' - conv(x1,x2))+...
    norm(squeeze(auto_der.derivative_Flambda(:,2,:)).' - 0)+...
    norm(squeeze(auto_der.derivative_Flambda(:,3,:)).' - 0);

max_error(2) = norm(auto_der.derivative_Fx_diagonal - [1,0,0;0,-4,0;0,0,1]);
max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix(1,1,:)).' - l1* conv(x2,one_vec)) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,1,:)).' - l1* conv(x1,one_vec) - conv(x3,one_vec)) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,1,:)).' - conv(x2,one_vec));

max_error(3) = max_error(3) + norm(squeeze(auto_der.derivative_Fx_toeplix(1,2,:)).' - conv(x3,one_vec) + conv(one_vec,one_vec)) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,2,:)).' - 0) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,2,:)).' - conv(x1,one_vec));

max_error(3) = max_error(3) + norm(squeeze(auto_der.derivative_Fx_toeplix(1,3,:)).' - conv(one_vec,one_vec)) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,3,:)).' - 0) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,3,:)).' - 0);


for i = 1:3
    for j = 1:3
        for k = 1:3
            if i == 3 && k == 3
                continue
            end
            if ~isempty(auto_der.derivative_Fx_delay{i,j}{k})
                error('wrong size delay derivative')
            end
        end
    end
end

max_error(4) = 0;
max_error(5) = 0;

if any(max_error > 10^-10)
    error('Test 2 delay derivative - 3D to 3D')
end
% string_vector_field = 'dot x1 + l1 * x1 * x2 + x2 * x3 \n -4 *dot x2 + x3 * x1 - x1 \n dot x3 + x1';

auto_mat_der = derivative_to_matrix(auto_der);

scal_coefs1 = conv(squeeze(scalar_equation.linear_coef{2}(:,1,:)).',one_vec);
scal_coefs2 = conv(squeeze(scalar_equation.linear_coef{2}(:,2,:)).',one_vec);
scal_coefs3 = conv(squeeze(scalar_equation.linear_coef{2}(:,3,:)).',one_vec);
scal_coef = [scal_coefs1, scal_coefs2, scal_coefs3];
% delay_mat = 0;
dot_der = kron([1,0,0; 0,-4,0;0,0,1], mat_bigK * 1i);
mat_conv = [ l1*conv_mat_big(conv(x2,one_vec)), conv_mat_big(conv(one_vec,l1*x1+x3)) , conv_mat_big(conv(x2,one_vec)) 
    conv_mat_big(conv(x3 - one_vec,one_vec)), conv_mat_big(conv(0*one_vec,one_vec)) , conv_mat_big(conv(x1,one_vec)) 
    conv_mat_big(conv(one_vec,one_vec)) , conv_mat_big(conv(0*x1,one_vec)) , conv_mat_big(conv(x2*0,one_vec)) ];
d_lambda = [conv(x1, x2), 0*conv(x1, x2), 0*conv(x1, x2)].';

hand_solution = [0, scal_coef;
    d_lambda, dot_der + mat_conv];

max_error = norm(hand_solution - auto_mat_der);

if max_error > 10^-10
    error('Test 2 delay derivative - 3D to 3D')
end



string_vector_field = 'dot x1 + l1 * x1 * x2 + x2^2 * x3 \n -4*dot x2 + x3 * x1^2 - x1 \n dot x3 + x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda(:,1,:)).' - conv(conv(x1,x2),one_vec))+...
    norm(squeeze(auto_der.derivative_Flambda(:,2,:)).' - 0)+...
    norm(squeeze(auto_der.derivative_Flambda(:,3,:)).' - 0);

max_error(2) = norm(auto_der.derivative_Fx_diagonal - [1,0,0;0,-4,0;0,0,1]);

max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix(1,1,:)).' - l1* conv(conv(x2,one_vec),one_vec)) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,1,:)).' - l1* conv(conv(x1,one_vec),one_vec) - 2*conv(conv(x3,x2),one_vec)) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,1,:)).' - conv(conv(x2,x2),one_vec));

max_error(3) = max_error(3) + norm(squeeze(auto_der.derivative_Fx_toeplix(1,2,:)).' - 2*conv(conv(x3,one_vec),x1) + conv(one_vec,conv(one_vec,one_vec))) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,2,:)).' - 0) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,2,:)).' - conv(x1,conv(x1,one_vec)));

max_error(3) = max_error(3) + norm(squeeze(auto_der.derivative_Fx_toeplix(1,3,:)).' - conv(one_vec,conv(one_vec,one_vec))) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,3,:)).' - 0) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,3,:)).' - 0);


for i = 1:3
    for j = 1:3
        for k = 1:3
            if i == 3 && k == 3
                continue
            end
            if ~isempty(auto_der.derivative_Fx_delay{i,j}{k})
                error('wrong size delay derivative')
            end
        end
    end
end

max_error(4) = 0;
max_error(5) = 0;

if max(max_error) > 10^-10
    error('Test 2b derivative delay - 3D to 3D')
end

% string_vector_field = 'dot x1 + l1 * x1 * x2 + x2^2 * x3 \n -4*dot x2 + x3 * x1^2 - x1 \n dot x3 + x1';
auto_mat_der = derivative_to_matrix(auto_der);
one_vec_big = conv(one_vec,one_vec);

scal_coefs1 = conv(squeeze(scalar_equation.linear_coef{2}(:,1,:)).',one_vec_big);
scal_coefs2 = conv(squeeze(scalar_equation.linear_coef{2}(:,2,:)).',one_vec_big);
scal_coefs3 = conv(squeeze(scalar_equation.linear_coef{2}(:,3,:)).',one_vec_big);
scal_coef = [scal_coefs1, scal_coefs2, scal_coefs3];
delay_mat = 0;
dot_der = kron([1,0,0; 0,-4,0;0,0,1], mat_verybigK * 1i);
mat_conv = [ l1*conv_mat_big(conv(x2,one_vec_big)), conv_mat_big(conv(conv(one_vec,l1*x1)+2*conv(x2,x3),one_vec)) , conv_mat_big(conv(conv(x2,x2),one_vec)) 
    conv_mat_big(conv(conv(2*x3,x1)-one_vec_big,one_vec)), conv_mat_big(conv(0*one_vec_big,one_vec)) , conv_mat_big(conv(conv(x1,x1),one_vec)) 
    conv_mat_big(conv(one_vec_big,one_vec)) , conv_mat_big(conv(0*one_vec_big,one_vec)) , conv_mat_big(conv(one_vec_big*0,one_vec)) ];
d_lambda = [conv(one_vec,conv(x1, x2)), 0*conv(one_vec_big, x2), 0*conv(one_vec_big, x2)].';

hand_solution = [0, scal_coef;
    d_lambda, dot_der + mat_conv];

max_error = norm(hand_solution - auto_mat_der);


if max(max_error) > 10^-10
    error('Test 2b derivative delay - 3D to 3D')
end


string_vector_field = 'dot x1 + l1 * Delay(x1,3.14) * x1 + x2^2 * x3 \n dot x2 + Delay(x2,3) * x1^2 - x1 \n l1 * dot x3 + x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda(:,1,:)).' - conv(conv(Delay(x1,3.14), x1),one_vec) - l1 * conv(1i *K*3.14.*Delay(x1,3.14),conv(x1,one_vec)))+...
    norm(squeeze(auto_der.derivative_Flambda(:,2,:)).' - conv(1i *K*3.*Delay(x2,3),conv(x1,x1)))+...
    norm(squeeze(auto_der.derivative_Flambda(:,3,:)).' - conv(conv(1i*K.*x3, one_vec),one_vec));

max_error(2) = norm(auto_der.derivative_Fx_diagonal - [1,0,0;0,1,0;0,0,l1]);

max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix(1,1,:)).' - l1* conv(conv(Delay(x1,3.14),one_vec),one_vec)) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,1,:)).' - 2 * conv(conv(x2,x3),one_vec)) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,1,:)).' - conv(conv(x2,x2),one_vec));

max_error(3) = max_error(3) + norm(squeeze(auto_der.derivative_Fx_toeplix(1,2,:)).' - 2*conv(conv(Delay(x2,3),x1),one_vec) + conv(one_vec, conv(one_vec, one_vec))) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,2,:)).' - 0) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,2,:)).' - 0);

max_error(3) = max_error(3) + norm(squeeze(auto_der.derivative_Fx_toeplix(1,3,:)).' - conv(one_vec,conv(one_vec,one_vec))) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,3,:)).' - 0) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,3,:)).' - 0);


for i = 1:3
    for j = 1:3
        for k = 1:3
            if i == 3 && k == 3
                continue
            end
            if all([i,j,k] == [1,1,2]) || all([i,j,k] == [2,2,2])
                if isempty(auto_der.derivative_Fx_delay{i,j}{k})
                    disp([i,j,k])
                    error('wrong size delay derivative')
                end
            else
                if ~isempty(auto_der.derivative_Fx_delay{i,j}{k})
                    disp([i,j,k])
                    error('wrong size delay derivative')
                end
            end
        end
    end
end
% term 1,1,2 : l1 * Delay(x1,3.14) * x1
% term 2,2,2 : Delay(x2,3) * x1^2
Delay = @(x, tau) exp(1i*l1*(-(length(x)-1)/2:(length(x)-1)/2)*tau).*x;
max_error(4) = norm(auto_der.derivative_Fx_delay{1,1}{2}.convolution.' - l1* conv(x1, conv(one_vec,one_vec))) + ...
    norm(auto_der.derivative_Fx_delay{2,2}{2}.convolution.' - conv(x1, conv(x1,one_vec)));

max_error(5) = norm(auto_der.derivative_Fx_delay{1,1}{2}.exp_coef - 1i*3.14*l1) + ...
    norm(auto_der.derivative_Fx_delay{2,2}{2}.exp_coef - 1i*3*l1);

if max(max_error) > 10^-10
    error('Test 3 delay derivative - 3D to 3D')
end

% string_vector_field = 'dot x1 + l1 * Delay(x1,3.14) * x1 + x2^2 * x3 \n dot x2 + Delay(x2,3) * x1^2 - x1 \n l1 * dot x3 + x1';
auto_mat_der = derivative_to_matrix(auto_der);
one_vec_big = conv(one_vec,one_vec);

scal_coefs1 = conv(squeeze(scalar_equation.linear_coef{2}(:,1,:)).',one_vec_big);
scal_coefs2 = conv(squeeze(scalar_equation.linear_coef{2}(:,2,:)).',one_vec_big);
scal_coefs3 = conv(squeeze(scalar_equation.linear_coef{2}(:,3,:)).',one_vec_big);
scal_coef = [scal_coefs1, scal_coefs2, scal_coefs3];
dot_der = kron([1,0,0; 0,1,0;0,0,l1], mat_verybigK * 1i);
mat_conv = [ l1*conv_mat_big(conv(Delay(x1,3.14),one_vec_big)), conv_mat_big(conv(2*conv(x2,x3),one_vec)) , conv_mat_big(conv(conv(x2,x2),one_vec)) 
    conv_mat_big(conv(conv(2*Delay(x2,3),x1)-one_vec_big,one_vec)), conv_mat_big(conv(0*one_vec_big,one_vec)) , 0*conv_mat_big(conv(conv(x1,x1),one_vec)) 
    conv_mat_big(conv(one_vec_big,one_vec)) , conv_mat_big(conv(0*one_vec_big,one_vec)) , conv_mat_big(conv(one_vec_big*0,one_vec)) ];

d_lambda = [conv(conv(Delay(x1,3.14), x1),one_vec) + l1 * conv(1i *K*3.14.*Delay(x1,3.14),conv(x1,one_vec)),...
    conv(1i *K*3.*Delay(x2,3),conv(x1,x1)), conv(conv(1i*K.*x3, one_vec),one_vec)].';

delay_mat11 = conv_mat_big(l1*conv(x1,one_vec_big)) *  expm(1i*3.14*l1 * mat_verybigK);
delay_mat22 = conv_mat_big(conv(x1,conv(x1,one_vec))) *  expm(1i*3*l1 * mat_verybigK);
delay_mat = [delay_mat11, 0*delay_mat11, 0*delay_mat11
    0*delay_mat22, delay_mat22, 0*delay_mat22
    0*delay_mat22, 0*delay_mat22, 0*delay_mat22];

hand_solution = [0, scal_coef;
    d_lambda, dot_der + delay_mat + mat_conv];

max_error = norm(hand_solution - auto_mat_der);


if max(max_error) > 10^-10
    error('Test 3 derivative delay - 3D to 3D')
end


string_vector_field = 'dot x1 + l1 * Delay(x1,3.14) * Delay(x2,4) + x2^2 * x3 \n dot x2 + Delay(x1,2) * Delay(x2,3)* x3 - x1 \n dot x3 + x1';

vector_field = from_string_to_polynomial_coef(string_vector_field,1,3);

full_zero_finding_problem = full_problem(scalar_equation,vector_field);

auto_der = derivative(full_zero_finding_problem, x);

max_error(1) = norm(squeeze(auto_der.derivative_Flambda(:,1,:)).' - conv(conv(Delay(x1,3.14), Delay(x2,4)),one_vec) ...
    - l1*conv(conv(1i *K*3.14.*Delay(x1,3.14), Delay(x2,4)),one_vec)...
    - l1*conv(conv(Delay(x1,3.14), 1i *K*4.*Delay(x2,4)),one_vec))+...
    norm(squeeze(auto_der.derivative_Flambda(:,2,:)).' ...
    - conv(conv(1i *K*2.*Delay(x1,2), Delay(x2,3)),x3)...
    - conv(conv(Delay(x1,2), 1i *K*3.*Delay(x2,3)),x3))+...
    norm(squeeze(auto_der.derivative_Flambda(:,3,:)).' - 0);

max_error(2) = norm(auto_der.derivative_Fx_diagonal - [1,0,0;0,1,0;0,0,1]);

max_error(3) = norm(squeeze(auto_der.derivative_Fx_toeplix(1,1,:)).' - 0) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,1,:)).' - 2 * conv(conv(x2,x3),one_vec)) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,1,:)).' - conv(conv(x2,x2),one_vec));

max_error(3) = max_error(3) + norm(squeeze(auto_der.derivative_Fx_toeplix(1,2,:)).' + conv(one_vec, conv(one_vec, one_vec))) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,2,:)).' - 0) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,2,:)).' - conv(conv(Delay(x2,3), Delay(x1,2)),one_vec));

max_error(3) = max_error(3) + norm(squeeze(auto_der.derivative_Fx_toeplix(1,3,:)).' - conv(one_vec,conv(one_vec,one_vec))) +...
    norm(squeeze(auto_der.derivative_Fx_toeplix(2,3,:)).' - 0) + ...
    norm(squeeze(auto_der.derivative_Fx_toeplix(3,3,:)).' - 0);


for i = 1:3
    for j = 1:3
        for k = 1:3
            if i == 3 && k == 3
                continue
            end
            if all([i,j,k] == [1,1,2]) || all([i,j,k] == [1,2,2]) || all([i,j,k] == [2,1,2]) || all([i,j,k] == [2,2,2])
                if isempty(auto_der.derivative_Fx_delay{i,j}{k})
                    disp([i,j,k])
                    error('wrong size delay derivative')
                end
            else
                if ~isempty(auto_der.derivative_Fx_delay{i,j}{k})
                    disp([i,j,k])
                    error('wrong size delay derivative')
                end
            end
        end
    end
end
% term 1,1,2 : l1 * Delay(x1,3.14) * Delay(x2,4)
% term 2,2,2 : Delay(x2,3) * Delay(x1,2) *x3
max_error(4) = norm(auto_der.derivative_Fx_delay{1,1}{2}.convolution.' - l1* conv(Delay(x2,4), conv(one_vec,one_vec))) + ...
    norm(auto_der.derivative_Fx_delay{1,2}{2}.convolution.' - l1 * conv(Delay(x1,3.14), conv(one_vec,one_vec))) + ...
    norm(auto_der.derivative_Fx_delay{2,1}{2}.convolution.' - conv(Delay(x2,3), conv(x3,one_vec))) + ...
    norm(auto_der.derivative_Fx_delay{2,2}{2}.convolution.' - conv(Delay(x1,2), conv(x3,one_vec)));


max_error(5) = norm(auto_der.derivative_Fx_delay{1,1}{2}.exp_coef - 1i*3.14*l1) + ...
    norm(auto_der.derivative_Fx_delay{1,2}{2}.exp_coef - 1i*4*l1) + ...
    norm(auto_der.derivative_Fx_delay{2,1}{2}.exp_coef - 1i*2*l1) + ...
    norm(auto_der.derivative_Fx_delay{2,2}{2}.exp_coef - 1i*3*l1);
if max(max_error) > 10^-10
    error('Test 4 derivative delay - 3D to 3D')
end

% string_vector_field = 'dot x1 + l1 * Delay(x1,3.14) * Delay(x2,4) + x2^2 * x3 \n dot x2 + Delay(x2,3) * Delay(x1,2) *x3 - x1 \n dot x3 + x1';
auto_mat_der = derivative_to_matrix(auto_der);
one_vec_big = conv(one_vec,one_vec);

scal_coefs1 = conv(squeeze(scalar_equation.linear_coef{2}(:,1,:)).',one_vec_big);
scal_coefs2 = conv(squeeze(scalar_equation.linear_coef{2}(:,2,:)).',one_vec_big);
scal_coefs3 = conv(squeeze(scalar_equation.linear_coef{2}(:,3,:)).',one_vec_big);
scal_coef = [scal_coefs1, scal_coefs2, scal_coefs3];
dot_der = kron([1,0,0; 0,1,0;0,0,1], mat_verybigK * 1i);
mat_conv = [ 0*conv_mat_big(conv(Delay(x1,3.14),one_vec_big)), conv_mat_big(conv(2*conv(x2,x3),one_vec)) , conv_mat_big(conv(conv(x2,x2),one_vec)) 
    conv_mat_big(conv(-one_vec_big,one_vec)), conv_mat_big(conv(0*one_vec_big,one_vec)) , conv_mat_big(conv(conv(Delay(x2,3),Delay(x1,2)),one_vec)) 
    conv_mat_big(conv(one_vec_big,one_vec)) , conv_mat_big(conv(0*one_vec_big,one_vec)) , conv_mat_big(conv(one_vec_big*0,one_vec)) ];
d_lambda = [ conv(conv(Delay(x1,3.14), Delay(x2,4)),one_vec) ...
    + l1*conv(conv(1i *K*3.14.*Delay(x1,3.14), Delay(x2,4)),one_vec)...
    + l1*conv(conv(Delay(x1,3.14), 1i *K*4.*Delay(x2,4)),one_vec),...
    + conv(conv(1i *K*2.*Delay(x1,2), Delay(x2,3)),x3)...
    + conv(conv(Delay(x1,2), 1i *K*3.*Delay(x2,3)),x3),...
    0*conv(conv(x1,x1),x1)].';


delay_mat11 = conv_mat_big(l1*conv(Delay(x2,4),one_vec_big)) *  expm(1i*3.14*l1 * mat_verybigK);
delay_mat12 = conv_mat_big(l1*conv(Delay(x1,3.14),one_vec_big)) *  expm(1i*4*l1 * mat_verybigK);
delay_mat21 = conv_mat_big(conv(conv(Delay(x2,3),x3),one_vec)) *  expm(1i*2*l1 * mat_verybigK);
delay_mat22 = conv_mat_big(conv(Delay(x1,2),conv(x3,one_vec))) *  expm(1i*3*l1 * mat_verybigK);
delay_mat = [delay_mat11, delay_mat12, 0*delay_mat11
    delay_mat21, delay_mat22, 0*delay_mat22
    0*delay_mat22, 0*delay_mat22, 0*delay_mat22];

hand_solution = [0, scal_coef;
    d_lambda, dot_der + delay_mat + mat_conv];

max_error = norm(hand_solution - auto_mat_der);

if max(max_error) > 10^-10
    error('Test 4 derivative delay - 3D to 3D')
end

disp('The construction of the derivative operator of a 3D vector field of a 3D vector with delay and convolutions is successful')

