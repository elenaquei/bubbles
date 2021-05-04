% hypercahotic system, figure 

global nu
global use_intlab 
global RAD_MAX
global talkative
use_intlab = 0;
nu = 1.1;
RAD_MAX = 10^-4;
%if isempty(talkative)
    talkative = 1;
%end

try
    sup(intval(1));
catch
    addpath(genpath('../'))
    startintlab
    first_run =1;
end


% b=c<sqrt(2e)
% d neq 0
% e>0
e = 2;
b=1;
c=b;
d=10;

f = @fn_hyper; % We choose a map.
df = @derivatives_hyper;
ddf = @second_der_hyper;
dddf = @third_der_hyper;
N = 4;
DIM= N;
X = [ -1.015513726331095
  1.464130774004092
  10.099019513592780
 -10.099019513592783
  -1.000000000000000
  30.610405390670862
  -0.607676658921156
   0.030310136744525
  -0.048061582141663
   3.964552188224825
   2.839481206874421
   0.062830225198066
   0.282991177919101
   0.788681629272493];
phi = [0.015584636790005
   0.609679882074390
   0.384647202347126
   0.867038269926880];

[x_star,lambda_star,eigenvec,eigenval, stability] = ...
    algebraic_hopf(f,df,ddf,dddf,N,X,phi);



string_hyper = '- dot x1 +a x1 -ax2-x2 x3+x4  \n - dot x2 -b x2 +x1 x3 \n - dot x3 -cx3+Dx1+x1x2 \n - dot x4 -ex1-ex2'; % fancy hyper
string_hyper_vars = strrep(string_hyper, 'b' , num2str(b)); % plugging in b
string_hyper_vars = strrep(string_hyper_vars, 'c' , num2str(c)); % plugging in c
string_hyper_vars = strrep(string_hyper_vars, 'D' , num2str(d)); % plugging in d
string_hyper_vars = strrep(string_hyper_vars, 'e' , num2str(e)); % plugging in e

string_hyper_cont = strrep(string_hyper_vars, 'a', 'l2'); % setting a as the second scalar variable

% some elements useful for the computation and the validation
n_nodes = 15; % number of Fourier nodes used: small, since near the Hopf bifurcation is a circle
n_iter = 1000; % number of iterations
step_size = 10^-4; % initial step size (then adapted along the validation
s = 'Hopf_hyper'; % where the solutions are stored

vectorfield = strrep(string_hyper_cont, 'l1' , '');
vectorfield = strrep(vectorfield, 'l2' , 'l1');
% string defining the vector field of the Hopf normal form 

f = from_string_to_polynomial_coef(vectorfield); 
% trasnformation into a vectorfield that can be used

eigenval = conj(mid(eigenval));
eigenvec = conj(mid(eigenvec));
sign_FLC = stability;

% starting the continuation
[s, last_sol, first_sol] = continuation_Hopf( lambda_star, x_star, f, n_nodes, n_iter, step_size, s, eigenvec, eigenval, sign_FLC, 2, 0);

load(s)

figure('DefaultAxesFontSize',14)
plot(x_n, 'b', 'LineWidth',2)
for j = 1:4
    subplot(1,4,j);
    hold on
end
plot(x_0, 'g', 'LineWidth',2)
max_t = max(2*pi * abs(x_0.scalar(1)),2*pi * abs(x_n.scalar(1)));
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


