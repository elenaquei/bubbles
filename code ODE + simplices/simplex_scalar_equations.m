function [list_of_nodes] = simplex_scalar_equations(node_number, ...
            bool_Hopf, list_of_nodes)
% function [simplex, list_of_nodes] = simplex_scalar_equations(simplex, ...
%           bool_Hopf, list_of_nodes)

for i = node_number
    node_i = list_of_nodes{i};
    if node_i.problem.scalar_equations.number_equations_lin == 4
        % checking if the equations have already been added
        continue
    end
    node_i = node_scalar_equation(node_i, bool_Hopf);
    list_of_nodes{i} = node_i;
end
end

function node = node_scalar_equation(node, bool_Hopf)

F_problem = node.problem;
x = node.solution;

if bool_Hopf
    % adds two new scalar equations
    F_problem = F_update_Hopf(F_problem,x);
else
    F_problem = fancy_scalar_condition(x, F_problem);
end

F_problem = continuation_equation_simplex(F_problem,x);

node.problem = F_problem;
end

