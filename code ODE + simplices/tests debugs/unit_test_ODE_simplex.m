% lorenz 84
% alpha = 0.25
% beta = 1
% G = 0.25
% delta = 1.04
% gamma = 0.987
% F = 2
% T parameter, starts at 0.0568
%
% dX = - Y^2 - Z^2 - alpha X + alpha F - gamma U^2
% dY = XY - beta XZ - Y+G
% dZ = beta XY + XZ -Z
% dU = delta U + gamma UX +T

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
    sup(intval(1));
catch
    startintlab;
end
%clearvars

%n_nodes= 60;
%DIM = 1;

% some elements useful for the computation and the validation
n_nodes = 7; % number of Fourier nodes used: small, since near the Hopf bifurcation is a circle
n_iter = 8;
bool_Hopf = 1;
step_size = 0.9*10^-2; % initial step size (then adapted along the validation
save_file = 'unit_test_ODE_simplex'; % where the solutions are stored

f = @fn_Lorenz84; % We choose a map.
df = @derivatives_Lorenz84;
ddf = @second_der_Lorenz84;
dddf = @third_der_Lorenz84;
N = 4;
DIM= N;

phi = [0.347693611712440   0.362607042331847   0.999178390839916   0.649549106979193]'; %lorenx84

X1 = [0.044235923390488
    0.481528303993606
    0.018000000000000
    0.664382383273932
    0.945151384904814
    0.677630170977385
    0.748019155686685
    0.656991171337332
    0.725337297889905
    0.354892602474774
    0.341136902930980
    0.318118801277541
    0.435591195800222
    0.972046133308465]; % gives mu = 0.010900890394664

F = 2;
alpha = 0.25;
beta =1 ;
G = 0.25;
delta = 1.04;
gamma = 0.987;
% T_null = 0.0568;

string_lorenz84 = '- dot x1 -x2^2 -x3^2 -alpha x1 + alpha*F x1^0- gamma x4^2 \n - dot x2 + x1 x2 - beta x1 x3 - x2 + G x1^0 \n - dot x3 + beta x1 x2 + x1 x3 - x3 \n - dot x4 - delta x4 + gamma x1 x4 + T x1^0'; % lorenz84
string_lorenz84_vars = strrep(string_lorenz84, 'alpha*F' , num2str(alpha*F));
string_lorenz84_vars = strrep(string_lorenz84_vars, 'alpha' , num2str(alpha)); % plugging in alpha
string_lorenz84_vars = strrep(string_lorenz84_vars, 'F' , num2str(F)); % plugging in F
string_lorenz84_vars = strrep(string_lorenz84_vars, 'G' , num2str(G)); % plugging in G
string_lorenz84_vars = strrep(string_lorenz84_vars, 'delta' , num2str(delta)); % plugging in delta
string_lorenz84_vars = strrep(string_lorenz84_vars, 'gamma' , num2str(gamma)); % plugging in gamma

string_lorenz84_cont1 = strrep(string_lorenz84_vars, 'T', 'l2'); % setting pho as the second scalar variable
string_lorenz84_cont2 = strrep(string_lorenz84_cont1, 'beta' , 'l3');

lor_rhs = @(x,T,beta)[(-x(2)^2-x(3)^2-alpha*x(1) + alpha*F- gamma *x(4)^2) ; (+ x(1) *x(2) - beta *x(1)* x(3) - x(2) + G);(+ beta *x(1) *x(2) + x(1)* x(3) - x(3)) ;(- delta *x(4) + gamma *x(1) *x(4) + T)];

vectorfield = strrep(string_lorenz84_cont1, 'l1' , '');%'-dot x1 - x2 + l1 x1 - x1 ^ 3 - x1  x2 ^ 2\n- dot x2 + x1 + l1 x2 - x1 ^ 2 x2 - x2 ^ 3';
vectorfield = strrep(vectorfield, 'l2' , 'l1');
vectorfield = strrep(vectorfield, 'beta' , num2str(beta));
vectorfield2 = strrep(string_lorenz84_cont2, 'l1' , '');%'-dot x1 - x2 + l1 x1 - x1 ^ 3 - x1  x2 ^ 2\n- dot x2 + x1 + l1 x2 - x1 ^ 2 x2 - x2 ^ 3';
vectorfield2 = strrep(vectorfield2, 'l2' , 'l1');
vectorfield2 = strrep(vectorfield2, 'l3' , 'l2');
% string defining the vector field of the Hopf normal form

f_lor = from_string_to_polynomial_coef(vectorfield);
% transformation into a vectorfield that can be used

[x0,lambda0,eigenvec,eigenval, stability] = ...
    algebraic_hopf(f,df,ddf,dddf,N,X1,phi);

use_intlab = 0;
sol_Xi = Hopf_system_setup(lambda0, x0, f_lor, n_nodes,...
    eigenvec, eigenval, stability, step_size);
% this gives us the wrong shape of things
sol = Xi_vector([sol_Xi.scalar(1),lambda0, beta, sol_Xi.scalar(3:end)], sol_Xi.vector);

polynomial = from_string_to_polynomial_coef(vectorfield2);
big_Hopf = Taylor_series_Hopf(polynomial,n_nodes);

load('stored_lorenz84.mat')
sol_N = list_of_nodes{5}.solution;
big_Hopf = list_of_nodes{5}.problem;

big_Hopf = F_update_Hopf(big_Hopf,sol_N);

use_intlab = 0;
bool_validated = 1;
save_file = continuation_simplex(sol_N, big_Hopf,...
    n_iter, step_size, save_file, bool_Hopf, bool_validated);

bool_validated = 0;
save_file = continuation_simplex(sol_N, big_Hopf,...
    n_iter, step_size, save_file, bool_Hopf, bool_validated);

load(save_file)
bool_validated = 1;
[list_of_simplices, index_non_validated, Interval, Z0_iter, ...
    Z1_iter, Z2_iter, Y_iter] = a_posteriori_validations(list_of_simplices,...
    list_of_nodes, [1], bool_Hopf);