% a straight forward example for anyone interested in starting using this
% library with continuation

% let's take a simple example of a polynomial ODE having periodic solution:
% \dot x = y
% \dot y = mu (1 - x^2 ) y - x
% with initial conditions (x(0),y(0)) = (2,0), and mu unknown
%
% We want to construct a numerical branch and validate it

% definition of some global variables - very important bit
%global first_run
global nu
global RAD_MAX
global Display
global use_intlab
global talkative
use_intlab = 0;
Display = 0;
RAD_MAX = 10^-5;
talkative = 2;

%if isempty(first_run)
try 
    intval(1);
catch
    addpath(genpath('../'))
    startintlab
    %first_run =1;
end

nu = 1.05; 
mu = 0.05;
n_iter = 500;
s='vdp_cont2';
% DATA OF THE PROBLEM
% number of nodes requested in Fourier
number_of_nodes = 16;

% rhs of the vector field and initial conditions
initial_coord = [2,0];
rhs= @(t,x) [x(2);(1-x(1).^2)*x(2)-x(1)];
% known approximation of the period of the solution 
approx_period = 6.8; % even a very rough approximation should do

% string defining the vector field, compatible with the requests of the
% function from_string_to_polynomial_coef 
string_vector_field = 'dot x1 - l1 x2 \n dot x2 -mu l1  x2 +mu l1  x1^2 x2 + l1 x1'; 

string_fixed = strrep(string_vector_field, 'mu', num2str(mu)); % for point wise first system, plugging in mu
% continuation vector field, where mu is the second scalar variable
string_cont = strrep(string_vector_field, 'mu', 'l2');


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
vector_field = from_string_to_polynomial_coef(string_fixed);

% construnction of the full zero-finding problem from the scalar equation
% and the vector field
F = full_problem(scalar_equation, vector_field);

% refinement of the numerical solution 
xXi = Newton_2(xXi,F,30,10^-10);

% defining the continuation problem in pho, starting at pho_null
sol = Xi_vector([xXi.scalar,mu], xXi.vector);
scal_eq = default_scalar_eq(sol);
polynomial = from_string_to_polynomial_coef(string_cont); 
F_not_square = full_problem(scal_eq, polynomial);


% calculation of the tangent vector
DF0 = derivative_to_matrix(derivative(F_not_square,sol,0));
x_dot_0 = kernel(DF0);
% symmetrisation (x_dot determined up to a complex value)
[~,index] = max(abs(real(x_dot_0(1:sol.size_scalar))));
angle = atan( imag(x_dot_0(index))/real(x_dot_0(index)));
x_dot_0 = exp( - 1i * angle) * x_dot_0; % bringing x_dot_0 to be symmetric (by multiplication with the appropriate complex rotation)
x_dot_0 = Xi_vec2vec(symmetrise(vec2Xi_vec(x_dot_0,sol)));

if x_dot_0(2)<0
    x_dot_0 = -x_dot_0;
end

% launch the validation
[s, x_n] = continuation ( sol, F_not_square, n_iter, 10^-4, x_dot_0,s, 10^-10);