% Hamiltonian
% a = 2
% b = 3
% c = 1
% d = -1
% mu parameter, around 0
% 
% dX = y
% dY = z
% dZ = u
% dU = az + b x + c x^2 + d x^3 + mu y

% DOESN'T WORK: seems to be a degenerate Hopf + saddle node validation
% also fails if we skip the non-degeneracy condition

global nu
global use_intlab 
global talkative 
global RAD_MAX
if isempty(talkative)
    talkative = 1;
end
use_intlab = 0;
nu = 1.1;
RAD_MAX = 10^-4;

try 
    intval(1);
catch
    addpath(genpath('../'));
    addpath(genpath('../../'))
    startintlab;
end
%clearvars

a = 2;
b = 3;
c = 1;
d = -1;

X = @(x) x(1); Y =@(x) x(2); Z =@(x) x(3); U = @(x)x(4);

f = @(x,mu) [Y(x)
    Z(x)
    U(x)
    a*Z(x) + b*X(x) + c *X(x)^2 + d* X( x)^3 + mu* Y(x)];

df = @derivatives_Hamiltonian;
ddf = @second_der_Hamiltonian;
dddf = @third_der_Hamiltonian;
N = 4;
DIM= N;
X = [    0.0000
   -1.0000
    0.0000
         0
   -0.0000
         0
    2.4614
   -0.9413
   -2.4614
    0.9413
   -0.9413
   -2.4614
    0.9413
    2.4614];
phi = [0.4086
    0.0058
    0.5442
    0.3603];

[x_star,lambda_star,eigenvec,eigenval, stability] = ...
    algebraic_hopf(f,df,ddf,dddf,N,X,phi,0);


% lambda_star = 0.001 * lambda_star;

string_Hamiltonian = '- dot x1 +x2 \n - dot x2 + x3 \n - dot x3 + x4 \n - dot x4 + a x3 + b x1 + c x1^2 + d  x1^3+  mu x2'; 
string_Hamiltonian_vars = strrep(string_Hamiltonian, 'a' , num2str(a));
string_Hamiltonian_vars = strrep(string_Hamiltonian_vars, 'b' , num2str(b)); % plugging in alpha
string_Hamiltonian_vars = strrep(string_Hamiltonian_vars, 'c' , num2str(c)); % plugging in F
string_Hamiltonian_vars = strrep(string_Hamiltonian_vars, '+ d ' , num2str(d)); % plugging in beta

string_Hamiltonian_cont = strrep(string_Hamiltonian_vars, 'mu', 'l2'); % setting pho as the second scalar variable


% some elements useful for the computation and the validation
n_nodes = 5; % number of Fourier nodes used: small, since near the Hopf bifurcation is a circle
n_iter = 50; % number of iterations
step_size = 0.5*10^-2; % initial step size (then adapted along the validation
s = 'Hopf_Hamilton'; % where the solutions are stored

vectorfield = strrep(string_Hamiltonian_cont, 'l1' , '');
vectorfield = strrep(vectorfield, 'l2' , 'l1');
% string defining the vector field of the Hopf normal form 

f_ham = from_string_to_polynomial_coef(vectorfield); 
% transformation into a vectorfield that can be used

% validated values converted to doubles 
eigenval = conj(mid(eigenval));
eigenvec = conj(mid(eigenvec));
%lambda_star = lambda_star;
%x_star = x_star;
%sign_FLC = sign(mid(dummy.l1));
sign_FLC = stability;

% starting the continuation

lambda0 = lambda_star;
x0 = x_star;
f = f_ham;
h = step_size;
index_saddle = 1;

[s, last_sol, first_sol] = continuation_Hopf( lambda_star, x_star, f, n_nodes, n_iter, step_size, s, eigenvec, eigenval, sign_FLC, index_saddle, 0);


%[s, x_n] = continuation ( sol_2, big_Hopf, n_iter, h, x_dot_0,s, 10^-6, bool_Hopf, bool_fancy_scalar, bool_saddle,num_index_sad);

load(s)
figure('DefaultAxesFontSize',14)
plot(norm_x(1,:),norm_x(3,:),'b.-','MarkerSize',20)
% yTick_float = -0.01:0.02:0.13;
% set(gca,'YTick',yTick_float)
% set(gca,'YTickLabel',sprintf('%1.2f|',yTick_float))
% xTick_float = 1+(0:0.02:0.1);
% set(gca,'XTick',xTick_float)
% set(gca,'XTickLabel',sprintf('%1.2f|',xTick_float))
% axis([0.99*xTick_float(1) 1.01*xTick_float(end) 0.99*yTick_float(1) 1.01*yTick_float(end)])

ylabel('$a$','Interpreter','Latex','FontSize',20)
xlabel('$\tau$','Interpreter','Latex','FontSize',20)


figure('DefaultAxesFontSize',14)
plot(first_sol, 'g', 'LineWidth',2)
for j = 1:4
    subplot(1,4,j);
    hold on
end
plot(last_sol, 'b', 'LineWidth',2)
max_t = max(2*pi * x_0.scalar(1),2*pi * x_n.scalar(1));
for j = 1:4
    subplot(1,4,j);
    hold on
    xlim([0,max_t])
    ylim auto
    lim_y = ylim;
    yTick_float = sort([linspace(lim_y(1),...
        lim_y(2),6),0]);
    set(gca,'YTick',yTick_float)
    set(gca,'YTickLabel',sprintf('%1.2f|',yTick_float))
    xlabel('t','Interpreter','Latex')
    ylabel(sprintf('$u_%i$',j),'Interpreter','Latex')
end
set(gcf,'position',[100,200,900,400])