% unit tests for scalar_eq with non_computable elements
% needs to be associated with the file "non_computable_func_unit_test.m"
global use_intlab
global nu
global talkative
global norm_weight
norm_weight = [1,1,1,1,1,1,1,1,1,0.8,1.2].';
talkative = 0;
nu = 1.001;
use_intlab = 0;
Rmax = 10^-2;
n_modes = 20; % TOO BIG AND IT'S a problem for Y, too small and it's a problem for Z1
n_iter = 600;
step_size = 10^-3;
plotting_instructions = 5;
bool_Hopf = 1; % the Hopf blow up is already taken into account
bool_validated = 1;

load('./Kevins_code/point_candidate.mat')
save_file = 'DDE_example';

% extract info from point_candidate
z0 = X_ref(8:8+2*N);
z1 = X_ref(9+2*N:9+4*N);
z2 = X_ref(10+4*N:end);

x0 = X_ref(1);
p = X_ref(2);
R0 = X_ref(3);
psi = X_ref(4);
a = X_ref(5);
eta1 = X_ref(6);
eta2 = X_ref(7);
mu = R0*exp(-p*x0);

xi = Xi_vector([psi, x0, a, mu, R0, p, eta1, eta2],[z0.';z1.';z2.']);
xi = reshape(xi, n_modes);
n_scal = xi.size_scalar;
n_vec = xi.size_vector;
n_nodes = xi.nodes;
xi = symmetrise(xi);

% create non-computable equations
n_non_comp_eqs = 3;
f_non_comp = @(x) noncomputable_eqs_for_DDE(x);

% not blowed up: \dot y + y - R0 h(Delay(y, tau)) * y * (1 - y)

% with polynomial embedding and blow up, with y = x + az
% F1 = - z0 + mu * x *( 1-x) * Delay( z1, tau) + mu * Delay(z2, tau) * ( -a * z0^2 + (1-2x)*z0)
% \dot z0 = F1
% \dot z1 = - p * z2 * F1 + eta1
% \dot z2 = - a p z2 * F1 + eta2
% scalar unknowns: psi, x, a, mu, p, eta1, eta2
% we know mu = R0 * exp(-p x)

minus_F1 = '+ z0 - mu x Delay(z1,tau) + mu x^2 Delay(z1,tau) + a mu Delay(z2,tau) z0^2 - mu Delay(z2,tau) z0  + 2 x mu Delay(z2,tau) z0';
minus_F2 = '- p   z2 z0 + p   mu x Delay(z1,tau) z2 - p   mu x^2 Delay(z1,tau) z2 - p   a mu Delay(z2,tau) z0^2 z2 + p   z2 mu Delay(z2,tau) z0  - 2 p   x mu Delay(z2,tau) z0 z2 + eta1';
minus_F3 = '- a p z2 z0 + a p mu x Delay(z1,tau) z2 - a p mu x^2 Delay(z1,tau) z2 - a ^2 p mu Delay(z2,tau) z0^2 z2 + a p z2 mu Delay(z2,tau) z0  - 2 a p x mu Delay(z2,tau) z0 z2 + eta2';

human_vector_field = strcat('psi dot z0', minus_F1, '\n psi dot z1', minus_F2, '\n psi dot z2', minus_F3);

string_vector_field = strrep(human_vector_field, 'psi' , 'l1');
string_vector_field = strrep(string_vector_field, 'x' , 'l2');
string_vector_field = strrep(string_vector_field, 'a ' , 'l3');
string_vector_field = strrep(string_vector_field, 'a^' , 'l3^');
string_vector_field = strrep(string_vector_field, 'mu' , 'l4');
string_vector_field = strrep(string_vector_field, 'R0' , 'l5');
string_vector_field = strrep(string_vector_field, 'p' , 'l6');
string_vector_field = strrep(string_vector_field, 'eta1' , 'l7');
string_vector_field = strrep(string_vector_field, 'eta2' , 'l8');
string_vector_field = strrep(string_vector_field, 'z0' , 'x1');
string_vector_field = strrep(string_vector_field, 'z1' , 'x2');
string_vector_field = strrep(string_vector_field, 'z2' , 'x3');
string_vector_field = strrep(string_vector_field, 'tau' , num2str(-tau));

vector_field = from_string_to_polynomial_coef(string_vector_field,n_scal,n_vec);

scalar_eq_pol = '-x + mu * x - mu x^2';
scalar_eq_pol = strrep(scalar_eq_pol, 'x' , 'l2');
scalar_eq_pol = strrep(scalar_eq_pol, 'mu' , 'l4');
polynomial = from_string_to_polynomial_coef(scalar_eq_pol,n_scal,n_vec);
lin_coefficients = cell(3,1);
lin_coefficients{2} = zeros(0,n_vec, 2*n_nodes+1);
scalar_equation = scalar_eq(0,1,n_scal, n_vec, lin_coefficients, polynomial, n_non_comp_eqs, f_non_comp);

scalar_equation = F_update_Hopf(scalar_equation,xi);

zero_finding_problem = full_problem(scalar_equation,vector_field);
full_zero_finding_problem = continuation_equation_simplex(zero_finding_problem, xi);
y = apply(full_zero_finding_problem, xi);
if norm(norm(y)) > 10^-6
    error('Initial approx too bad')
end

DF_old = derivative_to_matrix(derivative(full_zero_finding_problem, xi,0));
f = @(x)Xi_vec2vec(apply(full_zero_finding_problem,vec2Xi_vec(x,xi)));
%DF_num = numerical_der(f, Xi_vec2vec(xi));

xi = Newton_2(xi,full_zero_finding_problem);

xi_0 = xi;
xi_1 = xi;
xi_2 = xi;
if 3 == 1
    use_intlab = 1;
    % Y = Y_delay_simplex(full_zero_finding_problem,full_zero_finding_problem,full_zero_finding_problem, xi_0, xi_1, xi_2);
    
    % Z0 = Z0_delay_simplex(full_zero_finding_problem,full_zero_finding_problem,full_zero_finding_problem, xi_0, xi_1, xi_2);
    
    Z1 = Z1_delay_simplex(full_zero_finding_problem,full_zero_finding_problem,full_zero_finding_problem, xi_0, xi_1, xi_2);
    
    Z2 = Z2_delay_simplex(full_zero_finding_problem,full_zero_finding_problem,full_zero_finding_problem, xi_0, xi_1, xi_2, Rmax);
    use_intlab = 0;
end

first_node = node(1,xi,zero_finding_problem);
simplices_set_up(first_node, zero_finding_problem, step_size);

bool_validated = 0;
% only numerical run, then we can select the simplices we are interested in
% save_file = continuation_simplex(xi, zero_finding_problem,...
%     n_iter, step_size, save_file, bool_Hopf, bool_validated, plotting_instructions);
% 
% load(save_file)

smallest_amplitude = abs(xi.scalar(3));
while smallest_amplitude > 10^-4
    step_size_local = smallest_amplitude;
    n_iter_around = 20;
    disp(xi.scalar(3))
    
    save_file = continuation_simplex(xi, zero_finding_problem,...
        n_iter_around, step_size_local, save_file, bool_Hopf, bool_validated, plotting_instructions);
    load(save_file)
    amplitude_vec = 0 * [1:length(list_of_nodes)];
    for i=1:length(list_of_nodes)
        node_i = list_of_nodes{i};
        amplitude_vec(i) = node_i.solution.scalar(3);
    end
    [smallest_amplitude, best_node_index] = min(abs(amplitude_vec));
    best_node = list_of_nodes{best_node_index};
    xi = best_node.solution;
end


while smallest_amplitude > 3*10^-6
    step_size_local = 0.7*smallest_amplitude;
    n_iter_closing = 10;
    xi = best_node.solution;
    disp(xi.scalar(3))
    
    save_file = continuation_simplex(xi, zero_finding_problem,...
        n_iter_closing, step_size_local, save_file, bool_Hopf, bool_validated, plotting_instructions);
    load(save_file)
    amplitude_vec = 0 * [1:length(list_of_nodes)];
    for i=1:length(list_of_nodes)
        node_i = list_of_nodes{i};
        amplitude_vec(i) = node_i.solution.scalar(3);
    end
    [smallest_amplitude, best_node_index] = min(abs(amplitude_vec));
    best_node = list_of_nodes{best_node_index};
end


for i = 1:50
    step_size_local = 0.7*smallest_amplitude;
    n_iter_closing = 10;
    xi = best_node.solution;
    disp(xi.scalar(2))
    
    save_file = continuation_simplex(xi, zero_finding_problem,...
        n_iter_closing, 10^-2, save_file, bool_Hopf, bool_validated, plotting_instructions);
    load(save_file)
    mu_vec = 0 * [1:length(list_of_nodes)];
    for i=1:length(list_of_nodes)
        node_i = list_of_nodes{i};
        mu_vec(i) = node_i.solution.scalar(2);
    end
    [smallest_amplitude, best_node_index] = min(abs(mu_vec));
    best_node = list_of_nodes{best_node_index};
end


xi = best_node.solution;

save_file = continuation_simplex(xi, zero_finding_problem,...
    n_iter, step_size, save_file, bool_Hopf, bool_validated, plotting_instructions);
load(save_file)

xlabel('$x$','Interpreter','Latex', 'FontSize', 20);
ylabel('amplitude','Interpreter','Latex', 'FontSize', 20);
zlabel('$\mu$','Interpreter','Latex', 'FontSize', 20);

[list_of_simplices, index_non_validated, Interval, Z0_iter, ...
    Z1_iter, Z2_iter, Y_iter] = a_posteriori_validations(list_of_simplices,...
    list_of_nodes, [], bool_Hopf);


function [F, dxF, dxF_mat, dxxF_norm] = noncomputable_eqs_for_DDE(xi)
% function [F, dxF, dxF_mat, dxxF] = noncomputable_eqs_for_DDE(xi)
%
% xi can be a Xi_vector or the norm of a Xi_vector
global use_intlab

if isa(xi, 'Xi_vector')
    sum_x = zeros(xi.size_vector + xi.size_scalar,1);
    if isintval(xi) || use_intlab
        sum_x = intval(sum_x);
    end
    for j = 1: xi.size_scalar
        sum_x(j) = xi.scalar(j);
    end
    for j = 1:xi.size_vector
        sum_x(xi.size_scalar + j) = sum(xi.vector(j,:));
    end
else
    sum_x = xi;
    if isintval(xi) || use_intlab
        sum_x = intval(sum_x);
    end
end



% psi, x0, a, mu, R0, p, eta1, eta2],[z0;z1;z2]
psi = sum_x(1);
x = sum_x(2);
a = sum_x(3);
mu = sum_x(4);
R0 =sum_x(5);
p =sum_x(6);
z0 = sum_x(end-2);
z1 = sum_x(end-1);
z2 = sum_x(end);

g = @(u,epsilon) (epsilon~=0)* (exp(-epsilon*u)-1)/epsilon + (epsilon==0)*(-u);
d1g = @(u,epsilon) (epsilon~=0)*(-exp(-epsilon*u)) + (epsilon==0)*(-1);
d2g = @(u,epsilon) (epsilon~=0)*exp(-epsilon*u)*(-u*epsilon-1+exp(epsilon*u))/epsilon^2 + (epsilon==0)*0;
d1d1g = @(u, epsilon) (epsilon~=0)* epsilon * exp(-epsilon*u) + (epsilon==0)*0;
d1d2g = @(u, epsilon) (epsilon~=0)* u * exp(-epsilon*u) + (epsilon==0)*0;
d2d2g = @(u, epsilon) (epsilon~=0)* exp(-epsilon*u)/(epsilon^3) * ...
    (u^2*epsilon^2 + 2*u*epsilon - 2 * exp(u*epsilon) + 2) + (epsilon==0)*0;

abs_d1d2g = @(u, epsilon) (epsilon~=0)* u * exp(epsilon*u) + (epsilon==0)*0;

abs_d2d2g_plus = @(u, epsilon) exp(epsilon*u)/(epsilon^3) * ...
    (u^2*epsilon^2 + 2*u*epsilon - 2 * exp(u*epsilon) + 2);
abs_d2d2g_minus = @(u, epsilon) exp(epsilon*u)/(epsilon^3) * ...
    (u^2*epsilon^2 - 2*u*epsilon - 2 * exp(-u*epsilon) + 2);

abs_d2d2g = @(u, epsilon) (epsilon~=0)* max(abs_d2d2g_plus(u,epsilon),...
    abs_d2d2g_minus(u,epsilon)) + (epsilon==0)*0;

abs_d2d2g_plus = @(u, epsilon) (epsilon~=0)* exp(-epsilon*u)/(epsilon^3) * ...
    (u^2*epsilon^2 + 2*u*epsilon - 2 * exp(u*epsilon) + 2) + (epsilon==0)*0;
abs_d2d2g_minus = @(u, epsilon) (epsilon~=0)* exp(epsilon*u)/(epsilon^3) * ...
    (u^2*epsilon^2 - 2*u*epsilon - 2 * exp(-u*epsilon) + 2) + (epsilon==0)*0;
abs_d2d2g = @(u, epsilon) max( abs(abs_d2d2g_minus(u, epsilon)), abs(abs_d2d2g_plus(u, epsilon)));

abs_d1d1g = @(u, epsilon) (epsilon~=0)* epsilon * exp(epsilon*u) + (epsilon==0)*0;
abs_d1g =  @(u,epsilon) (epsilon~=0)*(exp(epsilon*u)) + (epsilon==0)*(1);

F = [mu - R0*exp(-p*x); -z1 + g(p*z0, a) ; -z2 + exp(-a * p * z0)];
if nargout == 1
    return
end
DpsiF = [0;0;0];
DxF = [ R0*p*exp(-p*x);0;0];
DaF = [ 0; d2g(p * z0, a); - p * z0 * exp(-a * p * z0)];
DmuF = [ 1; 0; 0];
DR0F = [-exp(-p*x); 0; 0];
DpF = [ R0 * x * exp(-p*x); d1g(p*z0, a)*z0; - a * z0 * exp(-a*p*z0)];
DetaF = [0;0;0];
Dz0F = [0; d1g(p*z0, a) * p ; -a * p * exp(-a*p*z0)];
Dz1F = [0;-1;0];
Dz2F = [0;0;-1];

dxF = [DpsiF, DxF, DaF, DmuF, DR0F, DpF, DetaF, DetaF, Dz0F, Dz1F, Dz2F];
if nargout == 2
    return
end
if isa(xi, 'Xivector')
    x_one_vec = 1 + 0*xi.vector(1,:);
    
    dxF_mat = [dxF(:,1:xi.size_scalar), kron(dxF(:,(xi.size_scalar+1):end), x_one_vec)];
else
    dxF_mat = NaN;
end
if nargout == 3
    return
end

DpsiDF = [0;0;0];
DxDF = [p*(R0 * p + 1 + R0 * x) * exp(p*x) + R0*exp(p*x); 0; 0];
DaDF = [ 0; abs_d1d2g(p*z0,a)*(z0+p) + abs_d2d2g(p*z0,a); (z0 + p*z0*(p*z0+a*z0+a*p))*exp(a*p*z0)];
DmuDF = [ 0; 0; 0];
DR0DF = [(p+x)*exp(-p*x); 0; 0];
DpDF = [ x*(R0*p + 1 + R0*x)*exp(p*x) + R0*exp(p*x);
    abs_d1d1g(p*z0,a)*(z0+p)*z0+abs_d1g(p*z0,a) + abs_d1d2g(p*z0,a)*z0;
    (z0+a)*exp(a*p*z0) + a*z0*(p*z0+a*z0+a*p)*exp(a*p*z0)];
DetaDF = [0;0;0];
Dz0DF = [0; abs_d1d1g(p*z0,a)*(z0+p)*p + abs_d1g(p*z0,a)+abs_d1d2g(p*z0,a)*p;
    (p+a)*exp(a*p*z0) + (p*z0+a*z0+a*p)*(a*p)*exp(a*p*z0)];
Dz1DF = [0;0;0];
Dz2DF = [0;0;0];

dxxF_norm = abs([DpsiDF, DxDF, DaDF, DmuDF, DR0DF, DpDF, DetaDF, DetaDF, Dz0DF, Dz1DF, Dz2DF]);

end
