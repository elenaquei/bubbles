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
n_nodes = 10;
step_size = 10^-4;
n_iter =10;
save_file = 'saved elements/FHN_simplex_validation'; % path where the validation will be saved

% starting parameters
%  (alpha,I,epsilon,gamma) = (0.1, 0.4, 0.2480, 1) 
alpha = 0.1;
I_fix = 0.4;
epsilon_fix = 0.2480;
% gamma = 1; % will not be coded for simplicity

% construct the numerical solution with forward integration from known
% initial conditions
init_coord  = [0.5667, 0.4864];
approx_period = 14.7;

% right hand side
f=@(t,x)[x(1)*(x(1)-alpha)*(1-x(1))-x(2)+I_fix;epsilon_fix*( x(1)-x(2) )];

% forward integration
[tout, yout] = ode23s(f,linspace(0,approx_period,n_nodes*20),init_coord);


% transformation to Xi_vector from time series
x_sol = time_series2Xi_vec(tout,yout,n_nodes);


% definition of the vector field in the form of a string

% ?- \dot x1 + x1^2 - x1 \alpha - x1^3- \alphau x1 - x2 + I\n - \dot x2 + \epsilon x1 - \epsilon \gamma x2
%
string_FHN = '- dot x1 + l1 x1^2 - alpha l1 x1 - l1 x1^3 + alpha l1 x1^2 - l1 x2 + I \n - dot x2 + epsilon l1 x1 - epsilon l1 x2'; % general FHN
% plug in all the variables
string_FHN_vars = strrep(string_FHN, 'alpha' , num2str(alpha)); 
% fixed point vector field
string_FHN_fix = strrep(string_FHN_vars, 'I' , num2str(I_fix)); 
string_FHN_fix = strrep(string_FHN_fix, 'epsilon' , num2str(epsilon_fix));
% setting I and epsilon as scalar variables
% stupidly, the change to Hopf vector field does not expect the period to
% be there already. Why? no idea...
string_FHN_vars = strrep(string_FHN_vars, 'l1' , ''); 
string_FHN_vars = strrep(string_FHN_vars, 'I' , 'l1'); 
string_FHN_vars = strrep(string_FHN_vars, 'epsilon' , 'l2');

x_sol.vector(:,1:n_nodes-4) = 0;
x_sol.vector(:,n_nodes+5:end) = 0;
x_sol = symmetrise(x_sol);

% constructing the ODEs systems
% fixed pho
% scal_eq = fancy_scalar_condition(x_sol);
scal_eq = default_scalar_eq(x_sol);
polynomial_fix = from_string_to_polynomial_coef(string_FHN_fix);
F_fix = full_problem(scal_eq, polynomial_fix);

% refining the approximation with Newton method
x_sol_N = Newton_2(x_sol,F_fix,30,10^-7);

polynomial = from_string_to_polynomial_coef(string_FHN_vars);
big_Hopf = Taylor_series_Hopf(polynomial,n_nodes);

% TO DO: build the new solution
% that means: change the solution we just found into a Hopf solution: 
% scalar = [omega, I, epsilon, a, x], vector = (y-x)/a
x = mean(yout,1);
amplitude = norm(yout - x);
moved_vector = x_sol_N.vector;
moved_vector(:,n_nodes+1) = moved_vector(:,n_nodes+1) - x.';
vector = moved_vector/amplitude;
sol = Xi_vector([x_sol_N.scalar,I_fix, epsilon_fix, amplitude, x], vector);

big_Hopf = F_update_Hopf(big_Hopf,sol);

% test_Hopf = big_Hopf;
% test_Hopf = continuation_equation_simplex(test_Hopf,sol);
% sol_N = Newton_2(sol,test_Hopf,30,10^-7);

% launch the validation
bool_Hopf = 1;
[list_of_simplices,list_of_nodes] = continuation_simplex( sol, big_Hopf,...
    n_iter, step_size, save_file, bool_Hopf);
