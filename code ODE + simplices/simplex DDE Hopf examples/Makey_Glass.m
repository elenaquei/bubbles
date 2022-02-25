% Makey Glass Delay Differential equation

global use_intlab
global nu
global talkative
global norm_weight
norm_weight = [1,1,1,1,1,1,1,1].';
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
load('./simplex DDE Hopf examples/candidate_mackey_blowup_v2.mat')
eta = 0;
xXi = time_series2Xi_vec(tt,yy.',n_modes);

xi = Xi_vector([period, x, a, tau, delta, eta],xXi.vector);

xi = reshape(xi, n_modes);
n_scal = xi.size_scalar;
n_vec = xi.size_vector;
n_nodes = xi.nodes;
xi = symmetrise(xi);

% create non-computable equations
n_non_comp_eqs = 1;
f_non_comp = @(x) noncomputable_eqs_for_MG(x);

s = cell(6,1);
s_coef = [-1, 1, -3, -6, -4, +1];
s{1} = ' \delta \tau z_0';
s{2} = ' \tau Delay(z_0,-1)';
s{3} = ' \tau x^4 Delay(z_0,-1) Delay(z_1,-1)';
s{4} = ' x^3 a \tau Delay(z_0,-1) Delay(z_0,-1)Delay(z_1,-1)';
s{5} = ' x^2 a^2 \tau Delay(z_0,-1) Delay(z_0,-1)Delay(z_0,-1) Delay(z_1,-1)';
s{6} = ' x a^3 \tau Delay(z_0,-1) Delay(z_0,-1)Delay(z_0,-1)Delay(z_0,-1)Delay(z_1,-1)';

s1 = cell(8,1);
s1_coef=[4, 4, 12, 12, 12, 12, 4, 4];
s1{1} = 'a z_1^2 x^3 ';
s1{2} = 'a z_1^2 x^7 ';
s1{3} = 'a^2 z_1^2 x^2 z_0 ';
s1{4} = 'a^2 z_1^2 x^6 z_0 ';
s1{5} = 'a^3 z_1^2 x z_0^2 ';
s1{6} = 'a^3 z_1^2 x^5 z_0^2'; 
s1{7} = 'a^4 z_1^2 z_0^3 ';
s1{8} = 'a^4 z_1^2 x^4 z_0^3';

F2 = '';
for i = 1:6
    for j = 1:8
        coef = s_coef(i)*s1_coef(j);
        if coef>0
            coef_str = append('+', int2str(coef));
        else
            coef_str = int2str(coef);
        end
        F2 = append(F2, coef_str, s{i}, s1{j});
    end
    F2 = append (F2,'\\');
end
F2 = append('dot z_1 ', F2,'-eta');


F1 = '\dot z_0+ \delta \tau z_0 +  a^3 \tau x Delay(z_0,-1)Delay(z_0,-1)Delay(z_0,-1)\\Delay(z_0,-1) Delay(z_1,-1) + 4 a^2 \tau x^2 Delay(z_0,-1)Delay(z_0,-1)Delay(z_0,-1) Delay(z_1,-1) \\+ 6 a \tau x^3 Delay(z_0,-1)Delay(z_0,-1) Delay(z_1,-1) \\+ 3 \tau x^4 Delay(z_0,-1)Delay(z_1,-1) - \tau Delay(z_0,-1) Delay(z_1,-1)';
% human readable
 fprintf( strrep(strrep(strrep(F1, '\\' , '\n'),'\d','d'),'\t','t')); disp(' ')
 fprintf( strrep(strrep(strrep(F2, '\\' , '\n'),'\d','d'),'\t','t')); disp(' ')
F1 = strrep(F1, '\' , '');
F2 = strrep(F2, '\' , '');

vector_field = strcat(F1, '\n', F2);

string_vector_field = strrep(vector_field, 'psi' , 'l1');
string_vector_field = strrep(string_vector_field, 'x' , 'l2');
string_vector_field = strrep(string_vector_field, 'delta' , 'l5');
string_vector_field = strrep(string_vector_field, 'tau' , 'l4');
string_vector_field = strrep(string_vector_field, 'a ' , 'l3');
string_vector_field = strrep(string_vector_field, 'a^' , 'l3^');
string_vector_field = strrep(string_vector_field, 'eta' , 'l6');
string_vector_field = strrep(string_vector_field, 'z_0' , 'x1');
string_vector_field = strrep(string_vector_field, 'z_1' , 'x2');


vector_field = from_string_to_polynomial_coef(string_vector_field,n_scal,n_vec);

scalar_eq_pol = 'x - delta x - delta x^5';
scalar_eq_pol = strrep(scalar_eq_pol, 'x' , 'l2');
scalar_eq_pol = strrep(scalar_eq_pol, 'delta' , 'l5');
polynomial = from_string_to_polynomial_coef(scalar_eq_pol,n_scal,n_vec);
lin_coefficients = cell(3,1);
lin_coefficients{2} = zeros(0,n_vec, 2*n_nodes+1);
scalar_equation = scalar_eq(0,1,n_scal, n_vec, lin_coefficients, polynomial, n_non_comp_eqs, f_non_comp);

scalar_equation = F_update_Hopf(scalar_equation,xi);

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
x = sum_x(2);
a = sum_x(3);

z0 = sum_x(end-1);
z1 = sum_x(end);

g = @(x,a,u) (1+(x+a*u)^4)*(1+x^4) ;
dxg = @(x,a,u) 4*(x+a*u)^3*(1+x^4) + (1+(x+a*u)^4)*4*x^3;
dag = @(x,a,u) 4*(x+a*u)^3*u*(1+x^4);
dug = @(x,a,u) 4*(x+a*u)^3*a*(1+x^4);

F = g(x,a,z0)*z1 - 1;

if nargout == 1
    return
end
DxF = dxg(x,a,z0)*z1;
DaF = dag(x,a,z0)*z1;
Dz0F = dug(x,a,z0)*z1;
Dz1F = g(x,a,z0);

dxF = [0, DxF, DaF, 0, 0, 0, Dz0F, Dz1F];
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
dxdxg = @(x,a,z0) 4 *(3 *(1 + x^4) *(x + a *z0)^2 + 8 *x^3* (x + a *z0)^3 + 3 *x^2* (1 + (x + a* z0)^4));
dadxg = @(x,a,z0) 24 *z0 *(x + a* z0) *(1 + 7* x^4 + 8* a *x^3* z0 + 2 *a^2* x^2 * z0^2);
dadag = @(x,a,z0) 12 *(1 + x^4) *z0^2 *(x + a* z0)^2;
dadz0g = @(x,a,z0) 4 *(1 + x^4)* (x + a *z0)^2 *(x + 4 *a* z0);
dz0dz0g = @(x,a,z0) 48* a^2* (1 + x^4) *(x + a* z0)^2;
dz0dxg = @(x,a,z0) 16* a* (x + a* z0)^2 *(3 + 7* x^4 + 4 *a* x^3 *z0);

DxDF = z1*(dxdxg(abs(x), abs(a), abs(z0)) + dadxg(abs(x), abs(a), abs(z0)) + dz0dxg(abs(x), abs(a), abs(z0)));
DaDF = z1*(dadxg(abs(x), abs(a), abs(z0)) + dadag(abs(x), abs(a), abs(z0)) + dadz0g(abs(x), abs(a), abs(z0)));
Dz0DF = z1*(dz0dxg(abs(x), abs(a), abs(z0)) + dadz0g(abs(x), abs(a), abs(z0)) + dz0dz0g(abs(x), abs(a), abs(z0)));
Dz1DF = dxg(abs(x), abs(a), abs(z0)) + dag(abs(x), abs(a), abs(z0)) + dz0g(abs(x), abs(a), abs(z0));

dxxF_norm = abs([0, DxDF, DaDF, 0, 0, 0, Dz0DF, Dz1DF]);

end

