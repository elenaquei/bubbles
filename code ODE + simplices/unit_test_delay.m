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
exact_sol = 1i*K.*x1 + l1 *  conv(exp(1i * l1 * K * 3.14) .* x1,x1,'same');

max_error = max(abs(auto_sol(3:end-1) - exact_sol));
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
exact_sol = 1i*K.*x1 + l1 *  conv(exp(1i * l1 * K * 3.14) .* x1,x2,'same');

max_error = max(abs(auto_sol(3:end-1) - exact_sol));
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
% STILL RETURNS THE IFFT - but I don't seem to be able to update this... whyyyyy??

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


