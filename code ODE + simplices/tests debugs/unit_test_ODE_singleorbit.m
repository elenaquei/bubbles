% a straight forward example for anyone interested in starting using this
% library

% let's take a simple example of a polynomial ODE having periodic solution:
% \dot x = y
% \dot y = mu (1 - x^2 ) y - x
% with initial conditions (x(0),y(0)) = (2,0), setting mu = 1
%
% We want to construct a numerical solution and validate it

% definition of some global variables - very important bit
%global first_run
global nu
global RAD_MAX
global Display
global use_intlab
global norm_weight
norm_weight = [];
use_intlab = 0;
Display = 0;
RAD_MAX = 10^-8;

%if isempty(first_run)
try 
    intval(1);
catch
    addpath(genpath('../'))
    startintlab
    %first_run =1;
end

nu = 1.05; 

% DATA OF THE PROBLEM
% number of nodes requested in Fourier
number_of_nodes = 60;

% rhs of the vector field and initial conditions
initial_coord = [2,0];
rhs= @(t,x) [x(2);(1-x(1).^2)*x(2)-x(1)];
% known approximation of the period of the solution 
approx_period = 6.8; % even a very rough approximation should do

% string defining the vector field, compatible with the requests of the
% function from_string_to_polynomial_coef 
string_vector_field = 'dot x1 - l1 x2 \n dot x2 - l1 x2 + l1 x1^2 x2 + l1 x1'; 


% SETTING UP THE DATA IN RIGHT FORMAT
% computation of the numerical solution with forward integration and built
% in ode solver
[tout, yout] = ode45(rhs,[0:0.001:approx_period], initial_coord);

% transformation to Xi_vector from time series
xXi = time_series2Xi_vec(tout,yout,number_of_nodes);

% construction of default scalar equation based on the Xi_vector just
% generated
scalar_equation = default_scalar_eq(xXi,2,0.4);

% construction of the vector field from the string
vector_field = from_string_to_polynomial_coef(string_vector_field);

% construnction of the full zero-finding problem from the scalar equation
% and the vector field
F = full_problem(scalar_equation, vector_field);

% refinement of the numerical solution 
xXi = Newton_2(xXi,F,30,10^-10);

% VALIDATION
success = validation_orbit(F, xXi);
