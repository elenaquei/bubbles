% smallest test for computation and validation of DDEs (including Hopf
% expansion)
global use_intlab
global nu
global talkative
global norm_weight
norm_weight = [1,1,1,1,1,1,1,1,1,0.8,1.2].';
talkative = 3;
nu = 1.001;
use_intlab = 0;
Rmax = 10^-2;
n_modes = 4; 
n_iter = 1;
step_size = 10^-10;
plotting_instructions = 1;
bool_Hopf = 1; % the Hopf blow up is already taken into account
bool_validated = 1;

load('./Kevins_code/point_candidate.mat')
save_file = 'tests debugs/unit_test_DDE';

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

first_node = node(1,xi,zero_finding_problem);
simplices_set_up(first_node, zero_finding_problem, step_size);

% run AND validate
save_file = continuation_simplex(xi, zero_finding_problem,...
    n_iter, step_size, save_file, bool_Hopf, bool_validated, plotting_instructions);

% run and later validate
step_size = 10^-8;
save_file = continuation_simplex(xi, zero_finding_problem,...
    n_iter, step_size, save_file, bool_Hopf, 0, plotting_instructions);

load(save_file)

[list_of_simplices, index_non_validated, Interval, Z0_iter, ...
    Z1_iter, Z2_iter, Y_iter] = a_posteriori_validations(list_of_simplices,...
    list_of_nodes, [1], bool_Hopf,save_file);



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


function plot(list_simplex, list_of_nodes, index_simplices, varargin)
if nargin <3 || isempty(index_simplices)
    index_simplices = 1:length(list_simplex);
end

for index = 1: length(index_simplices)
    i = index_simplices(index);
    if length(index_simplices)>10
        plot_simplex(list_simplex.simplex{i}, list_of_nodes, varargin{1},'LineStyle','none',varargin{2:end});
    else
        plot_simplex(list_simplex.simplex{i}, list_of_nodes, varargin{:});
    end
    hold on
end
alpha 0.5
set(gca,'FontSize',18)
hold off
xlabel('$x$','Interpreter','Latex', 'FontSize', 20);
ylabel('amplitude','Interpreter','Latex', 'FontSize', 20);
zlabel('$\mu$','Interpreter','Latex', 'FontSize', 20);
end

function plot_simplex(simplex, list_of_nodes, color, varargin)
% function plot_simplex(simplex, color)
if nargin<3 || isempty(color)
    if simplex.verified
        color = 'b';
    elseif simplex.frontal
        color = 'y';
    else
        color = 'r';
    end
end
p_coord = zeros(1,3);
R0_coord = zeros(1,3);
a_coord = zeros(1,3);
label = cell(3,1);
for i = 1:3
    p_coord(i) = list_of_nodes{simplex.nodes_number(i)}.solution.scalar(6);
    R0_coord(i) = list_of_nodes{simplex.nodes_number(i)}.solution.scalar(5);
    a_coord(i) = list_of_nodes{simplex.nodes_number(i)}.solution.scalar(3);
    label{i} = 'x'+string(simplex.nodes_number(i));
end
% label_simplex = 'S' +string(simplex.number);
% center_x = mean(p_coord);
% center_y = mean(R0_coord);
% center_z = mean(a_coord);
fill3(p_coord,R0_coord,a_coord,color, varargin{:})

%text(x_coord,y_coord,z_coord,label,'VerticalAlignment','bottom','HorizontalAlignment','right')
%text(center_x,center_y,center_z,label_simplex,'VerticalAlignment','bottom','HorizontalAlignment','right')
end
