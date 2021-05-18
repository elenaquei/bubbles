function [list_of_simplices,list_of_nodes] = continuation_simplex(x0, F,...
    n_iter, step_size, save_file, bool_Hopf)
% function list_of_simplices = continuation_simplex(x0, F, n_iter, h, save_file, bool_Hopf)

Interval = zeros(2,n_iter);
norm_x = zeros(x0.size_scalar+x0.size_vector,n_iter+6);
Z0_iter = zeros(x0.size_scalar+x0.size_vector,n_iter+6);
Z1_iter = zeros(x0.size_scalar+x0.size_vector,n_iter+6);
Z2_iter = zeros(x0.size_scalar+x0.size_vector,n_iter+6);
Y_iter = zeros(x0.size_scalar+x0.size_vector,n_iter+6);

%% set up the first 6 simplices
[list_of_simplices,list_of_nodes] = simplices_set_up(x0,F, step_size);

list_of_frontal_nodes = 1:6;

% validate all new simplices
for j = 1:6
    simplex = list_of_simplices.simplex{j};
    [flag,Imin,Imax,Yvector,Z0vector,Z1vector,Z2vector,...
        simplex, list_of_nodes] = ...
        radii_polynomials_simplex(simplex, list_of_nodes);
    if flag < 1
        error('Validation failed')
    end
    list_of_simplices.simplex{j} = simplex;
    %storage
    step_size(i) = h;
    norm_x(:,i)    = vert(norm(x1));
    norm_x(1:x0.size_scalar,i) = x1.scalar;% to avoid absolute values here
    Interval(:,i)  = [Imin,Imax]';
    Z0_iter(:,i)   = vert(Z0vector);
    Z1_iter(:,i)   = vert(Z1vector);
    Z2_iter(:,i)   = vert(Z2vector);
    Y_iter(:,i)    = vert(Yvector);
end



%% comtinuation proper
i = 0;
while i < n_iter
    % select node to grow
    node_number = select_growing_node(list_of_simplices,list_of_nodes,...
        list_of_frontal_nodes, F);
    
    % grow the node
    node = list_of_nodes{node_number};
    [index_new_simplices, list_of_nodes, list_of_simplices, ...
        list_of_new_frontal_nodes] = grow_simplex(node,...
        step_size, list_of_nodes, list_of_simplices, problem);
    
    if any(list_of_new_frontal_nodes<0)
        index = list_of_new_frontal_nodes<0;
        remove_node = abs(list_of_new_frontal_nodes(index));
        list_of_frontal_nodes(list_of_frontal_nodes==remove_node) = [];
        list_of_new_frontal_nodes(index)=[];
    end
    list_of_frontal_nodes = union(list_of_frontal_nodes,list_of_new_frontal_nodes);
    
    % add extra equations to the new simplices
    for k = index_new_simplices
        % TO DO !!!
        simplex = list_of_simplices.simplex{k};
        node_numbers = simplex.nodes_number;
        [simplex, list_of_nodes] = simplex_scalar_equations(node_numbers, ...
            bool_Hopf, list_of_nodes);
        list_of_simplices.simplex{k} = simplex;
        
    end
    
    % validate all new simplices
    for j = index_new_simplices
        simplex = list_of_simplices.simplex{j};
        [flag,Imin,Imax,Yvector,Z0vector,Z1vector,Z2vector,...
            simplex, list_of_nodes] = ...
            radii_polynomials_simplex(simplex, list_of_nodes);
        if flag < 1
            error('Validation failed')
        end
        list_of_simplices.simplex{j} = simplex;
        %storage
        step_size(j) = h;
        norm_x(:,j)    = vert(norm(x1));
        norm_x(1:x0.size_scalar,j) = x1.scalar;% to avoid absolute values here
        Interval(:,j)  = [Imin,Imax]';
        Z0_iter(:,j)   = vert(Z0vector);
        Z1_iter(:,j)   = vert(Z1vector);
        Z2_iter(:,j)   = vert(Z2vector);
        Y_iter(:,j)    = vert(Yvector);
    end
    i = i + length(index_new_simplices);
end

save(save_file,list_of_simplices,list_of_nodes,Interval,Z0_iter,Z1_iter,Z2_iter,Y_iter,norm_x,step_size);

end

function priority_node = select_growing_node(list_of_simplices, list_of_nodes, list_of_frontal_nodes, F)
gap_angle_vec = 0 * list_of_frontal_nodes;
for i=1:length(list_of_frontal_nodes)
    node_index = list_of_frontal_nodes(i);
    node = list_of_nodes{node_index};
    DF_x = derivative_to_matrix(derivative(F, node.solution,0));
    gap_angle_vec(i) = gap_angle(node, DF_x, list_of_simplices, list_of_nodes);
end

[~, priority_node_index] = min(abs(gap_angle_vec));
priority_node = list_of_frontal_nodes(priority_node_index);
min_gap = gap_angle_vec(priority_node_index);
end