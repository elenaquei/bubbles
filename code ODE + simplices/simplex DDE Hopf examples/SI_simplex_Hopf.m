% delay SI model 
%% global and importan variables 
global use_intlab
global nu
global talkative
global norm_weight
norm_weight = [1,1,1,1,1,1,1,1,1,1,1].';
talkative = 0;
nu = 1.001;
use_intlab = 0;
Rmax = 10^-2;
n_modes = 17; % too big it would be a problem for Y
step_size = 5*10^-6;
plotting_instructions = 50;
bool_Hopf = 1; % the Hopf blow up is already taken into account
bool_validated = 0;
n_iter = 400;

save_file = 'SI_model_data';

load('SI_start_at_Hopf.mat')
xi = reshape(xi, n_modes);
n_scal = xi.size_scalar;
n_vec = xi.size_vector;
n_nodes = xi.nodes;
xi = symmetrise(xi);
tau = 10;

%% building the vector field and the scalar equations

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

%% define the scalar equations (including the non-polynomial ones)
% create non-computable equations
n_non_comp_eqs = 3;
f_non_comp = @(x) noncomputable_eqs_for_DDE(x);
% polynomial equation
scalar_eq_pol = '-x + mu * x - mu x^2';
scalar_eq_pol = strrep(scalar_eq_pol, 'x' , 'l2');
scalar_eq_pol = strrep(scalar_eq_pol, 'mu' , 'l4');
polynomial = from_string_to_polynomial_coef(scalar_eq_pol,n_scal,n_vec);
% setting up scalar equations
lin_coefficients = cell(3,1);
lin_coefficients{2} = zeros(0,n_vec, 2*n_nodes+1);
scalar_equation = scalar_eq(0,1,n_scal, n_vec, lin_coefficients, polynomial, n_non_comp_eqs, f_non_comp);
% scalar equations for Hopf problems (include rescaling)
scalar_equation = F_update_Hopf(scalar_equation,xi);

zero_finding_problem = full_problem(scalar_equation,vector_field);
full_zero_finding_problem = continuation_equation_simplex(zero_finding_problem, xi);
y = apply(full_zero_finding_problem, xi);
if norm(norm(y)) > 10^-4
    error('Initial approx too bad')
end
%% starting the numerical search
xi = Newton_2(xi,full_zero_finding_problem);

first_node = node(1,xi,zero_finding_problem);
simplices_set_up(first_node, zero_finding_problem, step_size);


save_file = continuation_simplex(xi, zero_finding_problem,...
    n_iter, step_size, save_file, bool_Hopf, 0, plotting_instructions);
load(save_file)
%% validation - in patches to avoid memory overload
% plot_SI(list_of_simplices, list_of_nodes)
% n_simplices = length(list_of_simplices);
% for i = 1:n_simplices
%     [list_of_simplices, index_non_validated, Interval, Z0_iter, ...
%         Z1_iter, Z2_iter, Y_iter] = a_posteriori_validations(list_of_simplices,...
%         list_of_nodes, i, bool_Hopf,save_file);
%     % [list_of_simplices, index_non_validated, Interval, Z0_iter, ...
%     %    Z1_iter, Z2_iter, Y_iter] = a_posteriori_validations(list_of_simplices,...
%     %    list_of_nodes, [], bool_Hopf,save_file);
%     
%     save(save_file,'list_of_simplices','list_of_nodes','Interval','Z0_iter',...
%         'Z1_iter','Z2_iter','Y_iter','step_size','bool_Hopf', 'bool_validated',...
%         'list_of_frontal_nodes');
%     % guarantees you can stop the computations pretty much any time!
% end


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

g = @(u,epsilon) fun_g(u,epsilon) ;
d1g = @(u,epsilon) -exp(-epsilon*u);
d2g = @(u,epsilon) fun_d2g(u,epsilon);
%d1d1g = @(u, epsilon) epsilon * exp(-epsilon*u);
%d1d2g = @(u, epsilon) u * exp(-epsilon*u);
%d2d2g = @(u, epsilon) fun_d2d2g;

abs_d1d2g = @(u, epsilon) abs(u) * exp(epsilon*u);

abs_d2d2g = @(u, epsilon) abs(fun_d2d2g(u,epsilon));

abs_d1d1g = @(u, epsilon) abs(epsilon) * exp(epsilon*u);
abs_d1g =  @(u,epsilon) exp(epsilon*u);

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

function g = fun_g(u,epsilon)
if isintval(epsilon)
    if in(0,epsilon) || inf(abs(epsilon))<1E-5
        N=4;
        err = sup(abs(u*(epsilon*u)^(N+1)/factorial(N+2)*exp(abs(epsilon*u))));
        while sup(err)>1e-15
            N = N+5;
            err = sup(abs(u*(epsilon*u)^(N+1)/factorial(N+2)*exp(abs(epsilon*u))));
        end
        g_center = -u*( (intval(1)./factorial(1:N+1)) * (-epsilon*u).^(0:N).');
        if any(isnan(g_center)) % Fix elementwise power of complex intval bug.
            slow_intval_fix = intval(zeros(N+1,1));
            for j=1:N+1
                slow_intval_fix(j) = (-epsilon*u)^(j-1);
            end
            g_center = -u*( (intval(1)./factorial(1:N+1)) * slow_intval_fix);
        end
        g = g_center + midrad(0,err);
    else
        g = 1/epsilon*(exp(-epsilon*u)-1);
    end
else
    if abs(epsilon)<1E-5
        N = 4;
        err = abs(u*(epsilon*u)^(N+1)/factorial(N+2)*exp(abs(epsilon*u)));
        while sup(err)>1e-15
            N = N+5;
            err = abs(u*(epsilon*u)^(N+1)/factorial(N+2)*exp(abs(epsilon*u)));
        end
        g = -u*( (1./factorial(1:N+1)) * ((-epsilon*u).^(0:N).') );
    else
        g = 1/epsilon*(exp(-epsilon*u)-1);
    end
end
end

function d2g = fun_d2g(u,epsilon)
if isintval(epsilon)
    if in(0,epsilon) || inf(abs(epsilon))<1E-5
        N=4;
        err = sup(abs(u^2*(epsilon*u)^(N+1)/factorial(N+3)*exp(abs(epsilon*u))));
        while sup(err)>1e-15
            N = N+5;
            err = sup(abs(u^2*(epsilon*u)^(N+1)/factorial(N+3)*exp(abs(epsilon*u))));
        end
        d2g_center = u^2*( ((intval(1)./factorial(0:N))./(2:N+2)) * (-epsilon*u).^(0:N).');
        if any(isnan(d2g_center))   % Fix elementwise power of complex intval bug.
            slow_intval_fix = intval(zeros(N+1,1));
            for j=1:N+1
                slow_intval_fix(j) = (-epsilon*u)^(j-1);
            end
            d2g_center = u^2*( ((intval(1)./factorial(0:N))./(2:N+2)) * slow_intval_fix);
        end
        d2g = d2g_center + midrad(0,err);
    else
        d2g = (1/epsilon^2)*(1 - exp(-epsilon*u)*(epsilon*u+1));
    end
else
    if abs(epsilon)<1E-5
        N = 4;
        err = abs(u^2*(epsilon*u)^(N+1)/factorial(N+3)*exp(abs(epsilon*u)));
        while sup(err)>1e-15
            N = N+5;
            err = abs(u^2*(epsilon*u)^(N+1)/factorial(N+3)*exp(abs(epsilon*u)));
        end
        d2g = u^2*( ((1./factorial(0:N))./(2:N+2)) * (-epsilon*u).^(0:N).');
    else
        d2g = (1/epsilon^2)*(1 - exp(-epsilon*u)*(epsilon*u+1));
    end
end
end

function d2d2g = fun_d2d2g(u,epsilon)
if isintval(epsilon)
    if in(0,epsilon) || inf(abs(epsilon))<1E-5
        N=4;
        err = sup(abs(u^3*(epsilon*u)^(N+1)/factorial(N+4)*exp(abs(epsilon*u))));
        while sup(err)>1e-15
            N = N+5;
            err = sup(abs(u^2*(epsilon*u)^(N+1)/factorial(N+4)*exp(abs(epsilon*u))));
        end
        d2d2g_center = -u^3*( ((intval(1)./factorial(0:N))./(3:N+3)) * (-epsilon*u).^(0:N).');
        if any(isnan(d2d2g_center)) % Fix elementwise power of complex intval bug.
            slow_intval_fix = intval(zeros(N+1,1));
            for j=1:N+1
                slow_intval_fix(j) = (-epsilon*u)^(j-1);
            end
            d2d2g_center = -u^3*( ((intval(1)./factorial(0:N))./(3:N+3)) *  slow_intval_fix);
        end
        d2d2g = d2d2g_center + midrad(0,err);
    else
        d2d2g = (1/epsilon^3)*(-2 + exp(-epsilon*u)*(epsilon^2*u^2 + 2*epsilon*u + 2));
    end
else
    if abs(epsilon)<1E-5
        N = 4;
        err = abs(u^3*(epsilon*u)^(N+1)/factorial(N+4)*exp(abs(epsilon*u)));
        while sup(err)>1e-15
            N = N+5;
            err = abs(u^3*(epsilon*u)^(N+1)/factorial(N+4)*exp(abs(epsilon*u)));
        end
        d2d2g = -u^3*( ((1./factorial(0:N))./(3:N+3)) * (-epsilon*u).^(0:N).');
    else
        d2d2g = (1/epsilon^3)*(-2 + exp(-epsilon*u)*(epsilon^2*u^2 + 2*epsilon*u + 2));
    end
end
end