% validation of saddle nodes in a standard ODE system 
% EXAMPLE: Rychkov
%
% x_dot = y - x^5 + x^3 - lambda x
% y_dot = -x

% working as of 7th February 2019 


global nu
global use_intlab 
global talkative 
global RAD_MAX
global Display
talkative = 1;
Display = 0;
use_intlab = 0;
nu = 1.01;
RAD_MAX = 10^-2;


try 
    intval(1);
catch
    addpath(genpath('../'));
    addpath(genpath('../../'))
    startintlab;
end

n_nodes = 47;

step_size = 10^-3;
s = 'rychkov_test';
s_num = 'rychokov_num';
s_temp = 'temp';
min_res_N = 10^-9;
mu = 0.224;
f = @(t,x) [ (x(2)-(x(1).^5 - x(1)^3+mu *x(1))); -x(1)];

[~,y]=ode45(f,0:0.01:6.8,[1.3;-0.8]);
[t,y]=ode45(f,0:0.01:6.8,y(end,:));
x0 = time_series2Xi_vec(t,y,n_nodes);

%string_vf = 'dot x1 - l1  x2 + l1 x1^5 - 2.5 l1 x1^3 + l1 x1 \n dot x2 + l1 x1';
string_vf = 'dot x1 - l1  x2 + l1 x1^5 - l1 x1^3 + 0.224 l1 x1 \n dot x2 + l1 x1';
f = from_string_to_polynomial_coef(string_vf);

scalar_eqs_fixed = fancy_scalar_condition(x0);
% scalar_eqs_fixed = default_scalar_eq(x0,1);
% scalar_eqs_fixed.linear_coef{3} = -1;
F_fixed = full_problem(scalar_eqs_fixed,f);
x0_N = Newton_2(x0,F_fixed,[],min_res_N);

%string_vf = 'dot x1 - l1  x2 + l1 x1^5 - l2 l1 x1^3 + l1 x1 \n dot x2 + l1 x1';
string_vf = 'dot x1 - l1  x2 + l1 x1^5 - l1 x1^3 + l1 l2 x1 \n dot x2 + l1 x1';
f = from_string_to_polynomial_coef(string_vf);
x0_N.size_scalar = 2;
x0_N.scalar = [x0_N.scalar,mu];

scalarEqs =  fancy_scalar_condition(x0_N); %default_scalar_eq(x0_N,1);
F_not_square = full_problem(scalarEqs,f);

DF_not_square = derivative_to_matrix(derivative(F_not_square,x0_N,0));
x_dot = kernel(DF_not_square);
[~,index] = max(abs(real(x_dot(1:x0_N.size_scalar))));
angle = atan( imag(x_dot(index))/real(x_dot(index)));
x_dot = exp( - 1i * angle) * x_dot;
x_dot = Xi_vec2vec(symmetrise(vec2Xi_vec(x_dot,x0_N)));
if x_dot(2)>0
    x_dot = -x_dot;
end

n_iter = 17; % 16 would be much faster
h_num = 5*10^-3;
[s_num,x1] = continuation_numerical( x0_N, F_not_square, 13, h_num/3, x_dot,s_temp, min_res_N, 0 ,1); 
%load(s_num_temp)
%[s_num,x1] = continuation_numerical( x1, F_not_square, 8, h_num, -x_dot_1,s_num, min_res_N, 0 ,1); % 15 + 30 + 10
% getting closer to the saddle node
load(s_num)
no_fancy_condition = 0;
no_Hopf_system = 0;
check_for_saddle =1;
var_saddle = 2;
[s, x_n] = continuation ( x1, F_not_square, n_iter, step_size, -x_dot_1,s, ...
    min_res_N,no_Hopf_system,no_fancy_condition,check_for_saddle,var_saddle);% saddle at iteration 3


% h_num = 5*10^-3;
% [s_num,x1] = continuation_numerical( x0_N, F_not_square, 1, h_num/2, -x_dot,s_num, min_res_N, 0 ,1); %
% load(s_num)
% no_fancy_condition = 1;
% no_Hopf_system = 0;
% check_for_saddle =1;
% n_iter = 4000;
% scalarEqs =  fancy_scalar_condition(x1); 
% F_not_square = full_problem(scalarEqs,f);
% [s, x_n] = continuation ( x0_N, F_not_square, n_iter, h, x_dot_1,s, min_res_N,no_Hopf_system,no_fancy_condition,check_for_saddle);




load(s)



f1 = figure;
axes('FontSize',15)
plot(norm_x(2,:),max(norm_x(4,:),norm_x(3,:)),'g.-','MarkerSize',20)
hold on

x_saddle0 = saddle_x0_stored{1};
x_saddle1 = saddle_x1_stored{1};
norm_x0 = norm(x_saddle0);
norm_x1 = norm(x_saddle1);

plot(norm_x1(2,:),max(norm_x1(4,:),norm_x1(3,:)),'b.','MarkerSize',20)
plot(norm_x0(2,:),max(norm_x0(4,:),norm_x0(3,:)),'g.','MarkerSize',20)
xlabel('$\mu$','Interpreter','Latex','FontSize',18)
ylabel('max$\{\|v_1\|_\nu, \|v_2\|_\nu\}$','Interpreter','Latex','FontSize',18)
%plot(norm_x(2,end),norm_x(4,end),'k*','MarkerSize',20)
f2 = figure;
axes('FontSize',15)
plot2(x_n,'b','LineWidth',2)
hold on
plot2(x_0,'g','LineWidth',2)
plot2(saddle_x0_stored{1},'r','LineWidth',2)


