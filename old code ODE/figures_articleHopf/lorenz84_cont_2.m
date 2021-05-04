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
talkative = 0;
use_intlab = 1;
nu = 1.1;
RAD_MAX = 10^-4;

try 
    sup(intval(1));
catch
    %addpath(genpath('./'));
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
%X = [0.02 0.52   0.005   0.6   0.91   0.2   0.8976   0.7768   0.72   0.155   0.538   0.2578   0.47   0.51]'; %lorenx84
% gives mu = 0.056838929149399;
phi = [0.347693611712440   0.362607042331847   0.999178390839916   0.649549106979193]'; %lorenx84

X = [   0.044235923390488
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

[x_star,lambda_star,eigenvec,eigenval, stability] = ...
    algebraic_hopf(f,df,ddf,dddf,N,X,phi);



% X =[0.991380908070511
%    0.536571300769029
%    0.018000000000000
%    0.297367818824924
%    0.075579005137251
%    0.731181693699346
%    0.574541608287545
%    0.591657184180207
%    0.568456879079416
%    0.372376120925458
%    0.703809586107104
%    0.862916874311999
%    0.436074392210215
%    0.265907261089055]; % gives mu = -0.010901127516006

% X =[0.836509981142622
%    0.055244417920191
%    0.018000000000000
%    0.094291826193074
%    0.250532535955538
%    0.815534180223778
%    0.399267901357534
%    0.767400850441951
%    0.856330320249457
%    0.041620912834673
%    0.993634890078026
%    0.752255126627537
%    0.514729787861757
%    0.216805206774801]; % gives mu = -0.056838929149399;

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

lor_rhs = @(x1,x2,x3,x4,T)[(-x2^2-x3^2-alpha*x1 + alpha*F- gamma *x4^2) ; (+ x1 *x2 - beta *x1* x3 - x2 + G);(+ beta *x1 *x2 + x1* x3 - x3) ;(- delta *x4 + gamma *x1 *x4 + T)];

% some elements useful for the computation and the validation
n_nodes = 5; % number of Fourier nodes used: small, since near the Hopf bifurcation is a circle
n_iter_Hopf = 55;
step_size_scal = 10^-3; % initial step size (then adapted along the validation
s = 'Hopf_lorenz84'; % where the solutions are stored

vectorfield = strrep(string_lorenz84_cont, 'l1' , '');%'-dot x1 - x2 + l1 x1 - x1 ^ 3 - x1  x2 ^ 2\n- dot x2 + x1 + l1 x2 - x1 ^ 2 x2 - x2 ^ 3';
vectorfield = strrep(vectorfield, 'l2' , 'l1');
% string defining the vector field of the Hopf normal form 

f_lor = from_string_to_polynomial_coef(vectorfield); % trasnformation into a vectorfield that can be used

% definition of the solution
% x_star 
% lambda_star 

% validated values converted to doubles 
eigenval = conj(mid(eigenval));
eigenvec = conj(mid(eigenvec));
sign_FLC = stability;

% starting the continuation
index_saddle = 2;% only consider saddle with respect to the "real" 
                 % parameter (not the amplitude not the period)
bool_saddle = 0;
use_intlab = 1;
%[s, last_sol, first_sol] = continuation_Hopf( lambda_star, x_star, f_lor, ...
%    n_nodes, n_iter_Hopf, step_size, s, eigenvec, eigenval, sign_FLC, index_saddle, bool_saddle);

s_middle1 = 'lorenz84_close_to_Hopf1';
fprintf('Computing %s\n',s_middle1);

[s_middle1] = continuation_Hopf( lambda_star, x_star, f_lor, ...
    n_nodes, 10, step_size_scal/20, s_middle1, eigenvec, eigenval, sign_FLC, index_saddle, bool_saddle);

load(s_middle1)

bool_Hopf= 1; 
bool_fancy_scalar=1;

for iter1 = 2:60
    
    s_middle_part_iter = strcat('lorenz84_close_to_Hopf',sprintf('%i',iter1));
    
    fprintf('Computing %s\n',s_middle_part_iter);
    
    [s_middle_part_iter] = continuation ( x_n, F, 8, step_size_scal, x_dot_n, s_middle_part_iter,[],...
            bool_Hopf, bool_fancy_scalar);
    
    load(s_middle_part_iter)
end

fprintf('Testing if iteration %i could glue to the not-blowed-up system\n', iter1);

[y_n, y_n_dot, F_new] = hopf_gluing(s_middle_part_iter,f_lor,[],step_size(end)/150);

disp('completed first run')


%s_loc = s_middle_part_iter;
%use_intlab = 1;
%[y_n, y_n_dot, F_new] = hopf_gluing(s_loc,f_lor);

bool_Hopf= 0; 

s_second_part1 = 'Hopf_lorenz84_after_Hopf1';
lineLength = fprintf('Computing %s\n',s_second_part1);
step_size_scal = 2*10^-3;
[s_second_part1] = continuation ( y_n, F_new, 5, step_size_scal, y_n_dot,s_second_part1,...
        10^-6,  bool_Hopf, bool_fancy_scalar);
load(s_second_part1)


for iter2 = 2:20
    
    s_second_part_iter = strcat('Hopf_lorenz84_after_Hopf',sprintf('%i',iter2));
    
    lineLength = fprintf('Computing %s\n',s_second_part_iter);
    
    [s_second_part_iter] = continuation ( x_n, F, 5, min(step_size(end),step_size_scal), x_dot_n,s_second_part_iter,...
        10^-6, bool_Hopf, bool_fancy_scalar);
    
    load(s_second_part_iter)
end

disp('completed middle run')


fprintf('Testing if iteration %i could glue back to the blowed-up system', iter2);
[y_n, y_n_dot, F_new] = reverse_hopf_gluing(s_second_part_iter,s,f_lor,lor_rhs,[],step_size(end)/150);

%[y_n, y_n_dot, F_new] = reverse_hopf_gluing(s_second_part_iter,s,f_lor,lor_rhs);

% save('lorenz84_after_reverse_gluing')
% load('lorenz84_after_reverse_gluing')
clear saddle_x0_stored

%n_iter_last = 47;

bool_Hopf = 1;
bool_saddle = 0; % also validate the existence of a saddle node
bool_fancy_scalar = 1; % already included in the Hopf boolean
index_saddle = 2; 


iter3 = 1;
step_size_scal = 3*10^-3;
while iter3 < 150 && (~exist('saddle_x0_stored','var') || isempty(saddle_x0_stored) )
    
    s_third_part_iter = strcat('Hopf_lorenz84_to_second_Hopf',sprintf('%i',iter3));
    
    lineLength =  fprintf('Computing %s\n',s_third_part_iter);
    
    if iter3 == 1
        [s_third_part_iter] = continuation ( y_n, F_new, 5, step_size_scal, y_n_dot, s_third_part_iter,...
            10^-6,bool_Hopf, bool_fancy_scalar, bool_saddle, index_saddle);
    else
        [s_third_part_iter] = continuation ( x_n, F, 5, step_size_scal, x_dot_n, s_third_part_iter,...
            10^-6,bool_Hopf, bool_fancy_scalar, bool_saddle, index_saddle);
    end
    
    load(s_third_part_iter)
    
    iter3 = iter3 + 1;
end

iter3 = iter3 - 1;

if ~exist('saddle_x0_stored','var') || isempty(saddle_x0_stored)
    disp('DID NOT COMPLETE')
    disp('Did not reach the second Hopf bifurcation.')
else
    surface_plot_lorenz84(iter1, iter2, iter3)
end


return

load(s_second_part1)

figure('DefaultAxesFontSize',14)
plot(x_n, 'b', 'LineWidth',2)
for j = 1:4
    subplot(1,4,j);
    hold on
end
plot(x_0, 'g', 'LineWidth',2)
for j = 1:4
    subplot(1,4,j);
    hold on
    ylim auto
    lim_y = ylim;
    if ~any(linspace(lim_y(1),...
            lim_y(2),6)==0)
        yTick_float = sort([linspace(lim_y(1),...
            lim_y(2),6),0]);
    else
        yTick_float = linspace(lim_y(1),...
            lim_y(2),6);
    end
    set(gca,'YTick',yTick_float)
    set(gca,'YTickLabel',sprintf('%1.2f|',yTick_float))
    
end
set(gcf,'position',[100,200,900,400])
disp('initial solution in green, end solution in blu')

figure('DefaultAxesFontSize',14)
plot((norm_x(2,:)),norm_x(3:end,:).','.-','LineWidth',2)
xlabel('$\mu$','Interpreter','Latex','FontSize',20)
ylabel('$\|v_i\|_{\nu}$','Interpreter','Latex','FontSize',20)

