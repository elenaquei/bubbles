function [list_of_simplices,list_of_nodes] = continuation_simplex_Hopf(x0, F,...
    n_iter, step_size, save_file, Hopf_numerical)
% function list_of_simplices = continuation_simplex(x0, F, n_iter, h, save_file, bool_Hopf)

global use_intlab
global talkative
use_intlab = 0;

bool_Hopf = 1;

Interval = zeros(2,n_iter);
Z0_iter = zeros(x0.size_scalar+x0.size_vector,n_iter+6);
Z1_iter = zeros(x0.size_scalar+x0.size_vector,n_iter+6);
Z2_iter = zeros(x0.size_scalar+x0.size_vector,n_iter+6);
Y_iter = zeros(x0.size_scalar+x0.size_vector,n_iter+6);

%% set up the first 6 simplices
first_node = node(1,x0,F);
[list_of_simplices,list_of_nodes,list_of_frontal_nodes] = ...
    simplices_set_up(first_node, F, step_size);

plot(list_of_simplices);

% add the appropriate continuation equations
for k = 1:6
    simplex_k = list_of_simplices.simplex{k};
    node_numbers = simplex_k.nodes_number;
    for j = node_numbers
        n_near_nodes = neighboring_nodes(j, list_of_nodes, list_of_simplices);
        [list_of_nodes] = simplex_scalar_equations(j, ...
            bool_Hopf, list_of_nodes, n_near_nodes);
    end
end


% validate all new simplices
for j = 1:-6
    simplex_j = list_of_simplices.simplex{j};
    use_intlab = 1;
    %[flag,Imin,Imax,Yvector,Z0vector,Z1vector,Z2vector,...
    %    simplex_j, list_of_nodes] = ...
    %    radii_polynomials_simplex(simplex_j, list_of_nodes);
    use_intlab = 0;
    if flag < 1
        error('Validation failed')
    elseif talkative>0
        fprintf('The validation of the %i-th simplex succeeded\n', j)
        if talkative>1
            figure
            plot(list_of_simplices);
        end
    end
    list_of_simplices.simplex{j} = simplex_j;
    %storage
    %     Interval(:,j)  = [Imin,Imax]';
    %     Z0_iter(:,j)   = vert(Z0vector);
    %     Z1_iter(:,j)   = vert(Z1vector);
    Interval(:,j)  = [Imin,Imax]';
    Z0_iter(:,j)   = vert(Z0vector);
    Z1_iter(:,j)   = vert(Z1vector);
    Z2_iter(:,j)   = vert(Z2vector);
    Y_iter(:,j)    = vert(Yvector);
end



%% comtinuation proper
i = 0;
while i < n_iter
    % select node to grow
    node_number = select_growing_node(list_of_simplices,list_of_nodes,...
        list_of_frontal_nodes, F);
    fprintf('We grow from node %i\n',node_number)
    % grow the node
    node_i = list_of_nodes{node_number};
    % for debugging purposes
    list_of_simplices_old = list_of_simplices;
    list_of_nodes_old = list_of_nodes;
    
    [list_of_nodes,list_of_simplices, list_of_new_frontal_nodes, index_new_simplices]= ...
        grow_simplex(node_i, step_size, list_of_nodes, list_of_simplices, F);
    
    if mod(node_number,6)==0
    plot(list_of_simplices);
    pause
    end
    
    if any(list_of_new_frontal_nodes<0)
        index = list_of_new_frontal_nodes<0;
        remove_node = abs(list_of_new_frontal_nodes(index));
        for i = 1:length(remove_node)
            list_of_frontal_nodes(list_of_frontal_nodes==remove_node(i)) = [];
        end
        list_of_new_frontal_nodes(index)=[];
    end
    list_of_frontal_nodes = union(list_of_frontal_nodes,list_of_new_frontal_nodes);
    
    % add extra equations to the new simplices
    for index_k = 1:length(index_new_simplices)
        k = index_new_simplices(index_k);
        simplex = list_of_simplices.simplex{k};
        node_numbers = simplex.nodes_number;
        for jj = 1:3
            j = node_numbers(jj);
            n_near_nodes = neighboring_nodes(j, list_of_nodes, list_of_simplices);
            [list_of_nodes] = simplex_scalar_equations(node_numbers, ...
                bool_Hopf, list_of_nodes, n_near_nodes);
        end
        list_of_simplices.simplex{k} = simplex;
    end
    
    % validate all new simplices
    for index_j = 1:-length(index_new_simplices)
        j = index_new_simplices(index_j);
        simplex = list_of_simplices.simplex{j};
        [flag,Imin,Imax,Yvector,Z0vector,Z1vector,Z2vector,...
            simplex, list_of_nodes] = ...
            radii_polynomials_simplex(simplex, list_of_nodes);
        if flag < 1
            error('Validation failed')
        elseif talkative>0
            fprintf('The validation of the %i-th simplex succeeded\n', j)
            
        end
        list_of_simplices.simplex{j} = simplex;
        %storage
        Interval(:,j)  = [Imin,Imax]';
        Z0_iter(:,j)   = vert(Z0vector);
        Z1_iter(:,j)   = vert(Z1vector);
        Z2_iter(:,j)   = vert(Z2vector);
        Y_iter(:,j)    = vert(Yvector);
    end
    
    if talkative>1
        figure
        plot(list_of_simplices);
    end
    i = i + length(index_new_simplices);
end

save(save_file,list_of_simplices,list_of_nodes,Interval,Z0_iter,Z1_iter,Z2_iter,Y_iter,norm_x);

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
min_gap = gap_angle_vec(priority_node_index);
if abs(min_gap) > pi/5
    priority_node = min(list_of_frontal_nodes);
else
    priority_node = list_of_frontal_nodes(priority_node_index);
end

end

function list_near_nodes = neighboring_nodes(n_node, list_of_nodes, list_of_simplices)
patch = list_of_nodes{n_node}.patch;
list_near_nodes = [];
for index_j = 1: length(patch)
    j = patch(index_j);
    simplex_j = list_of_simplices.simplex{j};
    list_near_nodes = union(list_near_nodes, simplex_j.nodes_number);
end
list_near_nodes = setdiff( list_near_nodes, n_node);
if isempty(list_near_nodes)
    list_near_nodes = [];
end
end