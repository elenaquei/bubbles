global nu
global use_intlab
global talkative
global RAD_MAX
global Display
global norm_weight
norm_weight = [];
Display = 0;
talkative = 1;
use_intlab = 0;
RAD_MAX = 10^-2;


try
    intval(1);
catch
    startintlab
end

% problem dependent
nu = 1.05;
n_nodes = 7;
n_iter = 9000;
save_file = 'FHN_10K_25Jan'; % path where the validation will be saved

step_size = 0.01; % 0.7*10^-2
bool_Hopf = 1;
bool_validated = 0;
plotting_instructions = 50;

% starting parameters
%  (alpha,I,epsilon,gamma) = (0.1, 0.4, 0.2480, 1) 
alpha = 0.1;
I_fix = 0.4;
epsilon_fix = 0.2480;
% testing_epsilon_fix = .180;
% gamma = 1; % will not be coded for simplicity

% construct the numerical solution with forward integration from known
% initial conditions
init_coord  = [0.4491    0.4829]; % testing = [0.7099    0.4075];%
approx_period = 14.7;
% testing_approx_period = 16.74;

% right hand side
f=@(t,x)[x(1)^2 - alpha * x(1) - x(1)^3 + alpha * x(1)^2 - x(2) + I_fix ;
    epsilon_fix * x(1) - epsilon_fix * x(2)];

% definition of the vector field in the form of a string

% ?- \dot x1 + x1^2 - x1 \alpha - x1^3- \alphau x1 - x2 + I\n - \dot x2 + \epsilon x1 - \epsilon \gamma x2
%
string_FHN = '- dot x1 + l1 x1^2 - alpha l1 x1 - l1 x1^3 + alpha l1 x1^2 - l1 x2 + I l1 \n - dot x2 + epsilon l1 x1 - epsilon l1 x2'; % general FHN
% plug in all the variables
string_FHN_vars = strrep(string_FHN, 'alpha' , num2str(alpha)); 
% fixed point vector field
string_FHN_fix = strrep(string_FHN_vars, 'I' , num2str(I_fix)); 
string_FHN_fix = strrep(string_FHN_fix, 'epsilon' , num2str(epsilon_fix));
% setting I and epsilon as scalar variables
% stupidly, the change to Hopf vector field does not expect the period to
% be there already. Why? no idea...
string_FHN_vars = strrep(string_FHN_vars, 'l1' , ''); 
string_FHN_vars = strrep(string_FHN_vars, 'I' , 'l1'); 
string_FHN_vars = strrep(string_FHN_vars, 'epsilon' , 'l2');

% epsilon_vec = 0 * [1:length(list_of_nodes)];
% for i=1:length(list_of_nodes)
%     node = list_of_nodes{i};
%     epsilon_vec(i) = node.solution.scalar(3);
%     if node.solution.scalar(4)<0
%         epsilon_vec(i) =0;
%     end
% end
% [~, best_node_index] = max(epsilon_vec);
% best_node = list_of_nodes{best_node_index};
% 
% sol = best_node.solution;

% load initial solution
load('stored_FHN_start')
polynomial = from_string_to_polynomial_coef(string_FHN_vars);
big_Hopf = Taylor_series_Hopf(polynomial,n_nodes);

% build the new solution
% that means: change the solution we just found into a Hopf solution: 
% scalar = [omega, I, epsilon, a, x], vector = (y-x)/a

test_Hopf = big_Hopf;
test_Hopf = F_update_Hopf(test_Hopf,sol);
test_Hopf = continuation_equation_simplex(test_Hopf,sol);
sol_N = Newton_2(sol,test_Hopf,30,10^-7);

big_Hopf = F_update_Hopf(big_Hopf,sol_N);

% launch the validation
use_intlab = 0;
% save_file = continuation_simplex(sol_N, big_Hopf,...
%    n_iter, step_size, save_file, bool_Hopf, bool_validated, plotting_instructions);


% TOO MUCH memory usage, Matlab close to crashing: need to split the
% computation and the data into subsections
subsections = 200;
size_subsections = floor(n_iter/subsections);
n_simplices = n_iter;
last_index = 0;
for i = 1:subsections
    if i == subsections
        indices = last_index+1:n_simplices;
    else
        indices = last_index+(1:size_subsections);
    end
    last_index = last_index + size_subsections;
    load(save_file)
    [partial_list_of_simplices, partial_list_of_nodes] = ...
                subsample(list_of_simplices, indices, list_of_nodes);
    clear list_of_simplices list_of_nodes
    %plot(partial_list_of_simplices, partial_list_of_nodes)
    %hold
    save_file_iter = append('partial_FHN_simplex_validation',num2str(i));
    [partial_list_of_simplices_validated, index_non_validated, Interval, Z0_iter, ...
        Z1_iter, Z2_iter, Y_iter] = a_posteriori_validations(partial_list_of_simplices,...
        partial_list_of_nodes, [], bool_Hopf,[],save_file_iter);
end
