
global nu
global use_intlab
global talkative
global RAD_MAX
global Display
Display = 1;
talkative = 1;
use_intlab = 0;
RAD_MAX = 10^-2;

try
    intval(1);
catch
    startintlab
end

% problem dependent
nu = 1.05;
n_nodes = 50;
step_size = 10^-4;
n_iter =5;
save_file = 'saved elements/FHN_simplex_validation'; % path where the validation will be saved

% construct the numerical solution with forward integration from known
% initial conditions
init_coord  = [-1.2595115397689e+01  -1.6970525307084e+01   27];
approx_period = 2.3059072639398e+00;

% right hand side
f=@(t,x);

% forward integration
[tout, yout] = ode45(f,[0,approx_period],init_coord);

% transformation to Xi_vector from time series
xXi = time_series2Xi_vec(tout,yout,n_nodes);


% definition of the vector field in the form of a string

% ?- \dot x1 + x1^2 - x1 \alpha - x1^3- \alphau x1 - x2 + I\n - \dot x2 + \epsilon x1 - \epsilon \gamma x2
%
string_FHN = '- dot x1 + l1 x1^2 - l1 alpha x1 -l1  x1^3 - alpha l1 x1^2 - l1 x2 + I \n - dot x2 + epsilon l1 x1 - epsilon*gamma l1 x2'; % general FHN
% plug in all the variables
string_FHN_vars = strrep(string_FHN, 'epsilon*gamma' , num2str(epsilon*gamma));
string_FHN_vars = strrep(string_FHN_vars, 'alpha' , num2str(alpha)); 
string_FHN_vars = strrep(string_FHN_vars, 'I' , num2str(I)); 
% fixed point vector field
string_FHN_epsilon = strrep(string_FHN_vars, 'epsilon', num2str(epsilon)); 
% continuation vector field, where pho is the second scalar variable
string_FHN_cont = strrep(string_FHN_vars, 'gamma', num2str(gamma)); 
% setting pho as the second scalar variable

% constructing the ODEs systems
% fixed pho
scal_eq = default_scalar_eq(xXi);
polynomial_fix = from_string_to_polynomial_coef(string_FHN_epsilon);
F_fix = full_problem(scal_eq, polynomial_fix);

% refining the approximation with Newton method
xXi = Newton_2(xXi,F_fix,30,10^-7);

% defining the continuation problem in pho, starting at pho_null
sol = Xi_vector([xXi.scalar,pho_null], xXi.vector);
scal_eq = default_scalar_eq(sol);
polynomial = from_string_to_polynomial_coef(string_FHN_cont); 
F_not_square = full_problem(scal_eq, polynomial);


big_Hopf = Taylor_series_Hopf(polynomial,n_nodes);
% TO DO: build the new solution

% launch the validation
bool_Hopf = 1;
[list_of_simplices,list_of_nodes] = continuation_simplex ( sol, big_Hopf,...
    n_iter, step_size, save_file, bool_Hopf);
