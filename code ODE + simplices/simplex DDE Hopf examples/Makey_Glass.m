% Makey Glass Delay Differential equation

global use_intlab
global nu
global talkative
global norm_weight
norm_weight = [1,1,1,1,1].';
talkative = 0;
nu = 1.001;
use_intlab = 0;
Rmax = 10^-2;
n_modes = 17;
step_size = 5*10^-6;
plotting_instructions = 50;
n_iter = 4;

bool_Hopf = 1; 
bool_validated = 0;

save_file = 'Makey_Glass_data';

% extract info from point_candidate
load('./simplex DDE Hopf examples/candidate_mackey.mat')
eta = 0;
xXi = time_series2Xi_vec(tt,yy.',n_modes);

xi = Xi_vector([period, delta, tau],xXi.vector);

xi = reshape(xi, n_modes);
n_scal = xi.size_scalar;
n_vec = xi.size_vector;
n_nodes = xi.nodes;
xi = symmetrise(xi);

% create non-computable equations
n_non_comp_eqs = 1;
f_non_comp = @(x) noncomputable_eqs_for_MG(x);

vector_field = '- dot z_0 -delta tau z_0 + tau Delay(z_0,-1)Delay(z_1,-1) \n - dot z_1 + 4 delta tau z_0^4z_1^2 - 4 tau z_0^3z_1^2Delay(z_0,-1)Delay(z_1,-1)'; 

string_vector_field = strrep(vector_field, 'psi' , 'l1');
string_vector_field = strrep(string_vector_field, 'delta' , 'l2');
string_vector_field = strrep(string_vector_field, 'tau' , 'l3');
string_vector_field = strrep(string_vector_field, 'z_0' , 'x1');
string_vector_field = strrep(string_vector_field, 'z_1' , 'x2');

vector_field = from_string_to_polynomial_coef(string_vector_field,n_scal,n_vec);

lin_coefficients = cell(3,1);
lin_coefficients{2} = zeros(0,n_vec, 2*n_nodes+1);
scalar_equation = scalar_eq(0,0,n_scal, n_vec, lin_coefficients, polynomial_coefs(), n_non_comp_eqs, f_non_comp);

zero_finding_problem = full_problem(scalar_equation,vector_field);
full_zero_finding_problem = continuation_equation_simplex(zero_finding_problem, xi);
y = apply(full_zero_finding_problem, xi);
if norm(norm(y)) > 10^-4
    disp(norm(y))
    %error('Initial approx too bad')
end

DF_old = derivative_to_matrix(derivative(full_zero_finding_problem, xi,0));
f = @(x)Xi_vec2vec(apply(full_zero_finding_problem,vec2Xi_vec(x,xi)));

xi = Newton_2(xi,full_zero_finding_problem);
y = apply(full_zero_finding_problem, xi);
disp(norm(y))
disp(xi.scalar(6))

xi.scalar(6) = 0;
scalar_equation = F_update_Hopf(scalar_equation,xi);

zero_finding_problem = full_problem(scalar_equation,vector_field);
full_zero_finding_problem = continuation_equation_simplex(zero_finding_problem, xi);
xi = Newton_2(xi,full_zero_finding_problem);
y = apply(full_zero_finding_problem, xi);
disp(norm(y))
disp(xi.scalar(6))


first_node = node(1,xi,zero_finding_problem);
simplices_set_up(first_node, zero_finding_problem, step_size);

save_file = continuation_simplex(xi, zero_finding_problem,...
    n_iter, step_size, save_file, bool_Hopf, 0, plotting_instructions);
load(save_file)
[list_of_simplices, index_non_validated, Interval, Z0_iter, ...
Z1_iter, Z2_iter, Y_iter] = a_posteriori_validations(list_of_simplices,...
list_of_nodes, 1, bool_Hopf,[],save_file);

if bool_validated
    subsections = floor(n_iter/5);
    save_file_iter = 'MG/partial_MG_simplex_validation';
    validation_with_subpatches(save_file, save_file_iter, subsections)
end


function [F, dxF, dxF_mat, dxxF_norm] = noncomputable_eqs_for_MG(xi)
% function [F, dxF, dxF_mat, dxxF] = noncomputable_eqs_for_MG(xi)
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
z0 = sum_x(end-1);
z1 = sum_x(end);

g = @(z0, z1) z1 * (1+z0^4) - 1;

F = g(z0,z1);

if nargout == 1
    return
end
dz0g = @(z0,z1) z1*4*z0^3;
dz1g = @(z0,z1) 1+z0^4;

dxF = [0, 0, 0, dz0g(z0,z1), dz1g(z0,z1)];
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
% wolfram alpha results
dz0g = @(z0,z1) z1*4*z0^3;
dz1g = @(z0,z1) 1+z0^4;

dz0dz0g = @(z0,z1) z1*4*3*z0^2;
dz0dz1g = @(z0,z1) 4*z0^3;
dz1dz1g = @(z0,z1) 0;

Dz0DF = dz0dz0g(abs(z0),abs(z1)) + dz0dz1g(abs(z0),abs(z1));
Dz1DF = dz0dz1g(abs(z0),abs(z1)) + dz1dz1g(abs(z0),abs(z1));

dxxF_norm = abs([0, 0, 0, Dz0DF, Dz1DF]);

end

