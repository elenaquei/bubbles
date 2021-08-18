% lorenz 84
% alpha = 0.25
% beta = 1
% G = 0.25
% delta = 1.04
% gamma = 0.987
% F = 2
% T parameter, starts at 0.0568
%
% dX = - Y^2 - ^2 - alpha X + alpha F - gamma U^2
% dY = XY - beta XZ - Y+G
% dZ = beta XY + XZ -Z
% dU = delta U + gamma UX +T

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
n_iter = 35;
step_size = 0.5*10^-3; % initial step size 
save_file = 'Hopf_lorenz84_validated'; % where the solutions are stored
bool_validated = 1;

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

X2 =[0.991380908070511
    0.536571300769029
    0.018000000000000
    0.297367818824924
    0.075579005137251
    0.731181693699346
    0.574541608287545
    0.591657184180207
    0.568456879079416
    0.372376120925458
    0.703809586107104
    0.862916874311999
    0.436074392210215
    0.265907261089055]; % gives mu = -0.010901127516006


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
numerical_Hopfs = cell(2,1);
numerical_Hopfs{1}.x = x0;
numerical_Hopfs{1}.par_alpha = lambda0;
numerical_Hopfs{1}.par_beta = beta;
numerical_Hopfs{1}.eigenvec = eigenvec;
numerical_Hopfs{1}.eigenval = eigenval;

% numerically find the other Hopf solution! - I'm so glad I had the data!
% so so glad
[x1,lambda1,eigenvec1,eigenval1] = ...
    algebraic_hopf(f,df,ddf,dddf,N,X2,phi);
numerical_Hopfs{2}.x = x1;
numerical_Hopfs{2}.par_alpha = lambda1;
numerical_Hopfs{2}.par_beta = beta;
numerical_Hopfs{2}.eigenvec = eigenvec1;
numerical_Hopfs{2}.eigenval = eigenval1;

use_intlab = 0;
sol_Xi = Hopf_system_setup(lambda0, x0, f_lor, n_nodes,...
    eigenvec, eigenval, stability, step_size);
% this gives us the wrong shape of things
sol = Xi_vector([sol_Xi.scalar(1),lambda0, beta, sol_Xi.scalar(3:end)], sol_Xi.vector);

polynomial = from_string_to_polynomial_coef(vectorfield2);
big_Hopf = Taylor_series_Hopf(polynomial,n_nodes);

% big_Hopf.scalar_equations = big_Hopf_scalar_eqs(sol, numerical_Hopfs);
load('Hopf_lorenz84_numerical.mat')
sol_N = list_of_nodes{5}.solution;
big_Hopf_start = list_of_nodes{5}.problem;
% big_Hopf.scalar_equations = remove_lin_coef(big_Hopf.scalar_equations,6);
% big_Hopf.scalar_equations = remove_lin_coef(big_Hopf.scalar_equations,5);


% test_Hopf = big_Hopf;
% test_Hopf = F_update_Hopf(test_Hopf,sol);
% test_Hopf = continuation_equation_simplex(test_Hopf,sol);
% sol_N = Newton_2(sol,test_Hopf,30,10^-9);
% 
big_Hopf = F_update_Hopf(big_Hopf_start,sol_N);
% big_Hopf = continuation_equation_simplex(big_Hopf_start, sol_N, [1,2]);


bool_Hopf = 1;
use_intlab = 0;
save_file = continuation_simplex(sol_N, big_Hopf,...
    n_iter, step_size, save_file, bool_Hopf, bool_validated);


return

load(save_file)
plot(list_of_simplices,list_of_nodes)

new_name_file = 'new_Hopf_lorenz84_validated';
new_name_file = start_were_we_left_off(save_file, n_iter,new_name_file);

load(new_name_file)
plot(list_of_simplices,list_of_nodes)

new_name_file2 = 'new2_Hopf_lorenz84_validated';
new_name_file = start_were_we_left_off(new_name_file, n_iter, new_name_file2);

