% template_lor_cont
% 
% template to use correctly the continuation option of this library 

global nu
global use_intlab
global talkative
global RAD_MAX
global Display
global norm_weight
norm_weight = [];
Display = 0;
talkative = 1;
use_intlab = 0;
RAD_MAX = 10^-2;

try
    intval(1);
catch
    addpath(genpath('../'))
    startintlab
end

% problem dependent
nu = 1.05;
pho_null = 28;
n_nodes = 50;
step_size = 10^-4;
sigma = 10;
beta = 8/3;
pho = 28;
n_iter = 2;
s = 'tests debugs/unit_test_ODE_cont'; % path where the validation will be saved

% construct the numerical solution with forward integration from known
% initial conditions
init_coord  = [-1.2595115397689e+01  -1.6970525307084e+01   27];
approx_period = 2.3059072639398e+00;

% right hand side
f=@(t,x)[sigma*(x(2) - x(1)); x(1)*(pho-x(3)) - x(2); x(1)*x(2)- beta*x(3)];

% forward integration
[tout, yout] = ode45(f,[0,approx_period],init_coord);

% transformation to Xi_vector from time series
xXi = time_series2Xi_vec(tout,yout,n_nodes);


% definition of the vector field in the form of a string
string_lorenz = '- dot x1 + sigma l1 x2 - sigma l1 x1 \n - dot x2 + pho l1 x1 - l1 x1 x3 - l1 x2 \n - dot x3 + l1 x1 x2 - beta l1 x3'; % general lorenz
string_lorenz_vars = strrep(string_lorenz, 'sigma' , num2str(sigma)); % plugging in sigma
string_lorenz_vars = strrep(string_lorenz_vars, 'beta' , num2str(beta)); % plugging in beta
% fixed point vector field
string_lorenz_pho = strrep(string_lorenz_vars, 'pho', num2str(pho_null)); % for point wise first system, plugging in pho
% continuation vector field, where pho is the second scalar variable
string_lorenz_cont = strrep(string_lorenz_vars, 'pho', 'l2'); % setting pho as the second scalar variable

% constructing the ODEs systems
% fixed pho
scal_eq = default_scalar_eq(xXi);
polynomial_fix = from_string_to_polynomial_coef(string_lorenz_pho);
F_fix = full_problem(scal_eq, polynomial_fix);

% refining the approximation with Newton method
xXi = Newton_2(xXi,F_fix,30,10^-7);

% defining the continuation problem in pho, starting at pho_null
sol = Xi_vector([xXi.scalar,pho_null], xXi.vector);
scal_eq = default_scalar_eq(sol);
polynomial = from_string_to_polynomial_coef(string_lorenz_cont); 
F_not_square = full_problem(scal_eq, polynomial);


% launch the validation
[s, x_n] = continuation ( sol, F_not_square, n_iter, [], [],s, 10^-10);
