function save_file = continuation_simplex_core(save_file, list_of_simplices,list_of_nodes,...
    n_iter, step_size, bool_Hopf, bool_validated, list_of_frontal_nodes, ...
    Interval_old,Z0_iter_old, Z1_iter_old, Z2_iter_old, Y_iter_old)
global use_intlab
global talkative

if bool_validated
    Interval = zeros(2,n_iter);
    x0 = list_of_nodes{1}.solution;
    Z0_iter = zeros(x0.size_scalar+x0.size_vector,n_iter);
    Z1_iter = zeros(x0.size_scalar+x0.size_vector,n_iter);
    Z2_iter = zeros(x0.size_scalar+x0.size_vector,n_iter);
    Y_iter = zeros(x0.size_scalar+x0.size_vector,n_iter);
    if nargin > 7
        Interval = [ Interval_old, Interval];
        Z0_iter = [Z0_iter_old, Z0_iter];
        Z1_iter = [Z1_iter_old, Z1_iter];
        Z2_iter = [Z2_iter_old, Z2_iter];
        Y_iter = [Y_iter_old, Y_iter];
    end
end

F = list_of_nodes{end}.problem;

%% comtinuation proper
i = 0;
plot_index = [];
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
    
    if mod(node_number,30)==0
        hold on
        plot_index = union(plot_index,index_new_simplices);
        plot(list_of_simplices, list_of_nodes,plot_index);
        drawnow
        plot_index = [];
        if bool_validated
        save(save_file,'list_of_simplices','list_of_nodes','Interval',...
            'Z0_iter','Z1_iter','Z2_iter','Y_iter','step_size','bool_Hopf',...
            'bool_validated','list_of_frontal_nodes');
        else
            save(save_file,'list_of_simplices','list_of_nodes',...
                'step_size','bool_Hopf',...
                'bool_validated','list_of_frontal_nodes');
        end
        
        
    else
        plot_index = union(plot_index,index_new_simplices);
    end
    
    if any(list_of_new_frontal_nodes<0)
        index = list_of_new_frontal_nodes<0;
        remove_node = abs(list_of_new_frontal_nodes(index));
        for kj = 1:length(remove_node)
            list_of_frontal_nodes(list_of_frontal_nodes==remove_node(kj)) = [];
        end
        list_of_new_frontal_nodes(index)=[];
    end
    list_of_frontal_nodes = union(list_of_frontal_nodes,list_of_new_frontal_nodes);
    
    % validate all new simplices
    if bool_validated
        for index_j = 1:length(index_new_simplices)
            j = index_new_simplices(index_j);
            simplex_jj = list_of_simplices.simplex{j};
            use_intlab = 1;
            [flag,Imin,Imax,Yvector,Z0vector,Z1vector,Z2vector,...
                simplex_jj, list_of_nodes] = ...
                radii_polynomials_simplex(simplex_jj, list_of_nodes,bool_Hopf);
            use_intlab = 0;
            if flag < 1
                error('Validation failed')
            elseif talkative>0
                fprintf('The validation of the %i-th simplex succeeded\n', j)
            end
            simplex_jj.verified = 1;
            list_of_simplices.simplex{j} = simplex_jj;
            %storage
            Interval(:,j)  = [Imin,Imax]';
            Z0_iter(:,j)   = vert(Z0vector);
            Z1_iter(:,j)   = vert(Z1vector);
            Z2_iter(:,j)   = vert(Z2vector);
            Y_iter(:,j)    = vert(Yvector);
        end
    end
    
    %if talkative>1
    %    figure
    %    plot(list_of_simplices,list_of_nodes);
    %end
    i = i + length(index_new_simplices);
end

if bool_validated
save(save_file,'list_of_simplices','list_of_nodes','Interval','Z0_iter',...
    'Z1_iter','Z2_iter','Y_iter','step_size','bool_Hopf', 'bool_validated',...
    'list_of_frontal_nodes');
else
    save(save_file,'list_of_simplices','list_of_nodes',...
    'step_size','bool_Hopf', 'bool_validated',...
    'list_of_frontal_nodes');
end

end

function priority_node = select_growing_node(list_of_simplices, list_of_nodes, list_of_frontal_nodes, F)
gap_angle_vec = 0 * list_of_frontal_nodes;
epsilon_vec = 0 * list_of_frontal_nodes;
for i=1:length(list_of_frontal_nodes)
    node_index = list_of_frontal_nodes(i);
    node = list_of_nodes{node_index};
    DF_x = derivative_to_matrix(derivative(F, node.solution,0));
    gap_angle_vec(i) = gap_angle(node, DF_x, list_of_simplices, list_of_nodes);
    epsilon_vec(i) = node.solution.scalar(3);
end

[~, priority_node_index] = min(abs(gap_angle_vec));
min_gap = gap_angle_vec(priority_node_index);
if abs(min_gap) > pi/5
    priority_node = min(list_of_frontal_nodes);
    if priority_node > 100
       [~, min_eps_index] = max(abs(epsilon_vec));
       priority_node = list_of_frontal_nodes(min_eps_index);
    end
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
