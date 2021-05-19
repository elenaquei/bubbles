function [list_of_simplices,list_of_nodes, list_of_frontal_nodes] = ...
    simplices_set_up(starting_node,problem, step_size)
% function [list_of_simplices,list_of_nodes, list_of_frontal_nodes] = ...
%   simplices_set_up(x0,problem, step_size)
%
% INPUT
% node          instance of the node class - the node that needs "growing"
% step_size     expected length of the new edges
% list_of_nodes cell - all old nodes
% list_of_simplices     instance of the all_simplices class - list of all
%               the simplices (at leat all the once in the patch of node)
% problem       instance of the problem class - the ODE we are solving
%
% OUTPUT
% list_of_nodes     cell - new list of nodes
% list_of_simplices    instance of the all_simplices class - list of all
%                   the simplices (old + new)

list_of_nodes = cell(1,7);
list_of_nodes{1} = starting_node;

% wrapper
xc_h = step_size;
starting_solution = starting_node.solution;

div_tol = 1E-2; % Expansion/divergence criteria for step size reduction.
tol = 1E-12;    % Tolerance for convergence of Gauss-Newton.
F = @(X) Xi_vec2vec(apply(problem, vec2Xi_vec(X, starting_node)));
%F = @(X) X(1) - X(2)^2 - X(3)^2;
DF = @(X) derivative_to_matrix(derivative(problem,vec2Xi_vec(X, starting_solution),0));
% DF = @(X) [1, -2*X(2), -2*X(3)];


% Get projector
[Q,~] = qr(DF(Xi_vec2vec(starting_solution)).');
P = Q(:,end-1:end)*Q(:,end-1:end).'; % Projection to the tangent space

% Decide how many simplices are to be added to fill the gap.
n_new_simplex = 6;
gap = 2*pi;
y1_R2 = [0;1];
% Generate predictor "fan" in R2 coordinate system.
y_fan_R2 = zeros(2,n_new_simplex-1);
for j=1:n_new_simplex-1
    theta = gap*j/n_new_simplex;
    Rtheta = [cos(theta),sin(theta); -sin(theta),cos(theta)];
    
    y_fan_R2(:,j) = Rtheta*y1_R2;
    
end

R2_to_TM = null(DF(Xi_vec2vec(starting_solution)));
% Push the "fan" back into the tangent space T_{xc}M, and scale.
y_fan = R2_to_TM*y_fan_R2;
x_fan = zeros(dim,n_new_simplex-1);
for j=1:n_new_simplex-1
    y_fan(:,j) = y_fan(:,j)/norm(y_fan(:,j),2);
    x_fan(:,j) = starting_solution + y_fan(:,j)*xc_h;
end
% Refine predictors with Gauss-Newton
for j=1:n_new_simplex-1
    x_init = x_fan(:,j);
    h_local = xc_h;
    delta = inf;
    while delta>tol*(1+norm(x_init,2))
        if norm(x_init - x_fan(:,j))>div_tol   % Diverging.
            x_fan(:,j) = starting_solution + y_fan(:,j)/norm(y_fan(:,j),2)*h_local;
            x_init = x_fan(:,j);
        end
        x_fan(:,j) = GN(x_init,@(z)F(z),@(z)DF(z));
        delta = norm(x_fan(:,j)-x_init,2);
        x_init = x_fan(:,j);
    end
end

% transform x_new into a list of new simplices
list_new_nodes = zeros(1, n_new_simplex);
for i = 1:n_new_simplex
    number = i+1;
    solution = vec2Xi_vec(x_fan(:,i), starting_node);
    new_node = node(number, solution, problem);
    list_of_nodes{number} = new_node;
    list_new_nodes(i) = new_node.number;
end
frontal = 1;
verified = 0;
list_of_simplices = all_simplices();

for i = 1:n_new_simplex
    number1 = i+1;
    if i == n_new_simplex
        number2 = 2;
    else
        number2 = i;
    end
    nodes_number = [node.number, number1, number2];
    simplex_new = simplex(nodes_number, i, verified, frontal,list_of_nodes);
    list_of_simplices = append(list_of_simplices, simplex_new);
end
list_of_frontal_nodes = list_new_nodes;
list_new_nodes = union(list_new_nodes, out_nodes);
list_of_nodes = update_patches(list_of_nodes,list_of_simplices, list_new_nodes);


end


function x_new = GN(x,F,DF)
[Q,R] = qr(DF(x).');
[dim,~] = size(x);
R = R(1:dim-2,1:dim-2);
Q1 = Q(:,1:end-2);
z = R\F(x);
x_new = x - Q1*z;
end