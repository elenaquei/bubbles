global debug
debug = 0;
global nu
global use_intlab
global talkative
global RAD_MAX
global norm_weight
norm_weight = [];
nu = 1.1;
talkative = 2;
nu = 1.1;
RAD_MAX = 10^-4;

try
    infsup(1,2);
catch
    startintlab;
end

bool_validated = 1;
if ~bool_validated
    % some elements useful for the computation and the validation
    n_nodes = 3; % number of Fourier nodes used: small, since near the Hopf bifurcation is a circle
    n_iter = 70;
    % step size
    step_size = 0.5*10^-2; 
    save_file = 'Bautin_data'; % where the solutions are stored
else
    n_nodes = 3; 
    n_iter = 70;
    % step size
    step_size = 0.5*10^-2; 
    save_file = 'Bautin_data_validated'; % where the solutions are stored
end
save_file = 'Bautin_data';
if ~exist('Bautin', 'dir')
    mkdir('Bautin')
end
cd('Bautin')

N = 2;
DIM= N;

% l1 periodi
% l2 alpha
% l3 beta

string_Bautin = 'dot x - beta x + y - alpha x^3-alpha xy^2 + x^5 + 2x^3y^2 + xy^4 \n dot y - x - beta y - alpha yx^2 - alpha y^3 + yx^4 + 2y^3x^2 + y^5';

vectorfield = strrep(string_Bautin, 'l1' , '');
vectorfield = strrep(vectorfield, 'alpha' , 'l1');
vectorfield = strrep(vectorfield, 'beta' , 'l2');
vectorfield = strrep(vectorfield, 'x' , 'temp');
vectorfield = strrep(vectorfield, 'y' , 'x2');
vectorfield = strrep(vectorfield, 'temp' , 'x1');

f_Bautin = from_string_to_polynomial_coef(vectorfield);
% transformation into a vectorfield that can be used for computation and
% validation

use_intlab = 0;
x1 = zeros(1, 2 * n_nodes + 1);% sin
x1 ( n_nodes ) = 1/2*1i;
x1 (n_nodes + 2) = - 1/2*1i;
x2 = zeros(1, 2 * n_nodes + 1);
x2 ( n_nodes ) = 1/2;
x2 (n_nodes + 2) = 1/2;% cos
sol_N = Xi_vector([1/(2*pi), 0, 0, 0, 0, 0], [x1; x2]);


big_Hopf_start = Taylor_series_Hopf(f_Bautin,n_nodes);
big_Hopf = F_update_Hopf(big_Hopf_start,sol_N);

F_square = continuation_equation_simplex(big_Hopf, sol_N);
sol_N_better = Newton_2(sol_N,F_square);

bool_Hopf = 1;

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
