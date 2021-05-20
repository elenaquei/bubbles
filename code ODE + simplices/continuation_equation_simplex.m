function F_new = continuation_equation_simplex(F_new, x, ...
    list_of_nodes, n_near_nodes)

DF_x = derivative_to_matrix(derivative(F_new,x,0));
U = null(DF_x);
dim = length(x);

new_eq = F_new.scalar_equations.number_equations_lin+[1:2];

if nargin > 2
    other_basis = cell(length(n_near_nodes),0);
    counter = 1;
    for index = 1:length(n_near_nodes)
        j = n_near_nodes(index);
        node_nearby = list_of_nodes{j};
        local_scal_eqs = node_nearby.problem.scalar_equations;
        other_basis{counter} = zeros(2,dim);
        if local_scal_eqs.number_equations_lin > new_eq(1)
            other_basis{counter}(1,:) = extract_lin_coef(local_scal_eqs, new_eq(1));
            other_basis{counter}(2,:) = extract_lin_coef(local_scal_eqs, new_eq(2));
        end
        counter = counter + 1;
    end
    if ~isempty(other_basis)
        projection_U1 = zeros(length(other_basis),dim);
        projection_U2 = zeros(length(other_basis),dim);
        % project other_basis on U
        
        for i = 1: length(other_basis)
            projection_U1(i,:) = dot(other_basis{i}(1,:),U(:,1).')*U(:,1) ...
                + dot(other_basis{i}(1,:),U(:,2).')*U(:,2);
            projection_U2(i,:) = dot(other_basis{i}(2,:),U(:,1).')*U(:,1) ...
                + dot(other_basis{i}(2,:),U(:,2).')*U(:,2);
        end
        
        % average all projections
        average_proj1 = mean(projection_U1,1);
        average_proj2 = mean(projection_U2,1);
        U = [average_proj1;average_proj2].';
    end
end

lin_coef_vec1 = [ U(:,1).', - U(:,1).' * Xi_vec2vec(x)];

F_new.scalar_equations = change_lin_coef_vector(F_new.scalar_equations, ...
    lin_coef_vec1, new_eq(1)); 

lin_coef_vec2 = [ U(:,2).', - U(:,2).' * Xi_vec2vec(x)];
F_new.scalar_equations = change_lin_coef_vector(F_new.scalar_equations, ...
    lin_coef_vec2, new_eq(2)); 
end
