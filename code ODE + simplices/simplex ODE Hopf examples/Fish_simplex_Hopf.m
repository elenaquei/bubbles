% Church ODE: a pertubation ODE with two parameters, with epsilon far from
% 1 we might even have folds!
global debug
debug = 0;
global nu
global use_intlab
global talkative
global RAD_MAX
global norm_weight
norm_weight = [];
nu = 1.1;
talkative = 1;
nu = 1.1;
RAD_MAX = 10^-4;

try
    infsup(1,2);
catch
    startintlab;
end

bool_validated = 0;
if ~bool_validated
    % some elements useful for the computation and the validation
    n_nodes = 9; % number of Fourier nodes used: small, since near the Hopf bifurcation is a circle
    n_iter = 500;
    % step size
    step_size = 0.5*10^-2; 
    save_file = 'ChurchODE'; % where the solutions are stored
else
    n_nodes = 20; 
    n_iter = 5;
    % step size
    step_size = 10^-2; 
    save_file = 'ChurchODE_validated'; % where the solutions are stored
end
save_file = 'ChurchODE_test';

N = 3;
DIM= N;

epsilon = 0.8;
% l1 periodi
% l2 alpha
% l3 beta

string_Church = 'dot x - beta x + y +x^3+xy^2+xz^2+alpha^2x \n dot y - x - beta y + yx^2+y^3+epsilon yz^2+alpha^2y\n dot z + z^5 - 3z^3 + 0.01 z - 0.1 alpha';
string_Church_weird = 'dot x - beta x + y +x^3+xy^2+xz^2+alpha^2x \n dot y - x - beta y + yx^2+y^3+epsilon yz^2+alpha^2y\n dot z + z^5 - 3z^3 + 0.01 z - 0.1 alpha-0.01x^2-0.01y^2';

string_Church_vars = strrep(string_Church_weird, 'epsilon' , num2str(epsilon));

vectorfield = strrep(string_Church_vars, 'l1' , '');
vectorfield = strrep(vectorfield, 'alpha' , 'l1');
vectorfield = strrep(vectorfield, 'beta' , 'l2');
vectorfield = strrep(vectorfield, 'x' , 'temp');
vectorfield = strrep(vectorfield, 'z' , 'x3');
vectorfield = strrep(vectorfield, 'y' , 'x2');
vectorfield = strrep(vectorfield, 'temp' , 'x1');
% string defining the vector field of the Hopf normal form

f_Church = from_string_to_polynomial_coef(vectorfield);
% transformation into a vectorfield that can be used for computation and
% validation

use_intlab = 0;
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

% save_file = continuation_simplex(sol_N_better, big_Hopf,...
%     n_iter, step_size, save_file, bool_Hopf, bool_validated);
% new_name_file = 'ChurchODE_2_validated';
% new_name_file = start_were_we_left_off(save_file, n_iter,new_name_file);

save_file = continuation_simplex(sol_N_better, big_Hopf,...
    3*n_iter, step_size, save_file, bool_Hopf, 0);
xlabel('$\lambda_2$','Interpreter','Latex', 'FontSize', 20);
ylabel('amplitude','Interpreter','Latex', 'FontSize', 20);
zlabel('$\lambda_3$','Interpreter','Latex', 'FontSize', 20);

if bool_validated
    
    load(save_file)
    [list_of_simplices, index_non_validated, Interval, Z0_iter, ...
        Z1_iter, Z2_iter, Y_iter] =  a_posteriori_validations(list_of_simplices,...
        list_of_nodes, [], bool_Hopf, 1, save_file)
    
    save(save_file,'list_of_simplices','list_of_nodes','Interval','Z0_iter',...
        'Z1_iter','Z2_iter','Y_iter','step_size','bool_Hopf', 'bool_validated',...
        'list_of_frontal_nodes');
end
