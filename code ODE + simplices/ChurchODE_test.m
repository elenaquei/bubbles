% Church ODE: a pertubation ODE with two parameters, with epsilon far from
% 1 we might even have folds!
global debug
debug = 0;
global nu
global use_intlab
global talkative
global RAD_MAX
nu = 1.1;
talkative = 1;
nu = 1.1;
RAD_MAX = 10^-4;

try
    infsup(1,2);
catch
    startintlab;
end

% some elements useful for the computation and the validation
n_nodes = 9; % number of Fourier nodes used: small, since near the Hopf bifurcation is a circle
n_iter = 1000;
step_size = 5*10^-3; % initial step size 
save_file = 'ChurchODE'; % where the solutions are stored
bool_validated = 0;

N = 3;
DIM= N;

epsilon = 1;
% l1 periodi
% l2 alpha
% l3 beta

string_Church = '-dot x1 + l1 l3 x1 - l1 x2 - l1 x1 x1^2+ l1 x1 x2^2+ l1 x1 x3+l1 l2^2 x1\n -dot x2 + l1 x1 + l3 x2 -l1 x2 x1^2 -l1 x2 x2^2-epsilon l1 x2  x3-l1  l2^2 x2\n - dot x3 -l1 x3^5 + 3 l1 x3^3 - l1 x3 + l2'; % ChurchODE
string_Church_vars = strrep(string_Church, 'epsilon' , num2str(epsilon));

vectorfield = strrep(string_Church_vars, 'l1' , '');
vectorfield = strrep(vectorfield, 'l2' , 'l1');
vectorfield = strrep(vectorfield, 'l3' , 'l2');
% string defining the vector field of the Hopf normal form

f_Church = from_string_to_polynomial_coef(vectorfield);
% transformation into a vectorfield that can be used for computation and
% validation

use_intlab = 0;
n_nodes = 5;
x1 = zeros(1, 2 * n_nodes + 1);% sin
x1 ( n_nodes ) = 1/2*1i;
x1 (n_nodes + 2) = - 1/2*1i;
x2 = zeros(1, 2 * n_nodes + 1);
x2 ( n_nodes ) = 1/2;
x2 (n_nodes + 2) = 1/2;% cos
x3 = 0*x1;
sol_N = Xi_vector([1/(2*pi), 0, 0, 0, 0, 0, 0], [x1; x2; x3]);


big_Hopf_start = Taylor_series_Hopf(f_Church,n_nodes);
big_Hopf = F_update_Hopf(big_Hopf_start,sol_N);

F_square = continuation_equation_simplex(big_Hopf, sol_N);
sol_N_better = Newton_2(sol_N,F_square);

bool_Hopf = 1;

save_file = continuation_simplex(sol_N_better, big_Hopf,...
    n_iter, step_size, save_file, bool_Hopf, bool_validated);
