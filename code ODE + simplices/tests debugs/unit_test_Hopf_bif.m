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
global norm_weight
norm_weight = [];
nu = 1.1;
if isempty(talkative)
    talkative = 1;
end
use_intlab = 0;
nu = 1.1;
RAD_MAX = 10^-4;

try 
    intval(1);
catch
    addpath(genpath('./'));
    %addpath(genpath('../../'))
    startintlab;
end
%clearvars

%n_nodes= 60;
%DIM = 1;

f = @fn_Lorenz84; % We choose a map.
df = @derivatives_Lorenz84;
ddf = @second_der_Lorenz84;
dddf = @third_der_Lorenz84;
N = 4;
DIM= N;
X = [0.02 0.52   0.005   0.6   0.91   0.2   0.8976   0.7768   0.72   0.155   0.538   0.2578   0.47   0.51]'; %lorenx84
% gives mu = 0.056838929149399;
phi = [0.347693611712440   0.362607042331847   0.999178390839916   0.649549106979193]'; %lorenx84

[x_star,lambda_star,eigenvec,eigenval, stability] = ...
    algebraic_hopf(f,df,ddf,dddf,N,X,phi);

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
string_lorenz84_vars = strrep(string_lorenz84_vars, 'beta' , num2str(beta)); % plugging in beta
string_lorenz84_vars = strrep(string_lorenz84_vars, 'G' , num2str(G)); % plugging in G
string_lorenz84_vars = strrep(string_lorenz84_vars, 'delta' , num2str(delta)); % plugging in delta
string_lorenz84_vars = strrep(string_lorenz84_vars, 'gamma' , num2str(gamma)); % plugging in gamma

% string_lorenz84_pho = strrep(string_lorenz84_vars, 'T', num2str(T_null)); % for point wise first system, plugging in T
string_lorenz84_cont = strrep(string_lorenz84_vars, 'T', 'l2'); % setting pho as the second scalar variable


% some elements useful for the computation and the validation
n_nodes = 5; % number of Fourier nodes used: small, since near the Hopf bifurcation is a circle
n_iter = 3; % number of iterations
step_size = 10^-3; % initial step size (then adapted along the validation
s = 'tests debugs/unit_test_Hopf_bif'; % where the solutions are stored

vectorfield = strrep(string_lorenz84_cont, 'l1' , '');%'-dot x1 - x2 + l1 x1 - x1 ^ 3 - x1  x2 ^ 2\n- dot x2 + x1 + l1 x2 - x1 ^ 2 x2 - x2 ^ 3';
vectorfield = strrep(vectorfield, 'l2' , 'l1');
% string defining the vector field of the Hopf normal form 

f_lor = from_string_to_polynomial_coef(vectorfield); % trasnformation into a vectorfield that can be used

% definition of the solution
% x_star 
% lambda_star 

%dummy = load('hopf_in_lorenz84');

% validated values converted to doubles 
eigenval = conj(mid(eigenval));
eigenvec = conj(mid(eigenvec));
%lambda_star = lambda_star;
%x_star = x_star;
%sign_FLC = sign(mid(dummy.l1));
sign_FLC = stability;
%eigenvec = eigenvec / norm(eigenvec);
bool_saddle = 0;
index_saddle = 2;

% starting the continuation
[s, last_sol, first_sol] = continuation_Hopf( lambda_star, x_star, f_lor, ...
    n_nodes, n_iter, step_size, s, eigenvec, eigenval, sign_FLC,[],bool_saddle);
