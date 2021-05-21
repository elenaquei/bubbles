function [list_of_simplices,list_of_nodes] = continuation_simplex_Hopf(x0, F,...
    n_iter, step_size, save_file, rhs, phi, first_Hopf_data)
% function list_of_simplices = continuation_simplex(x0, F, n_iter, h, save_file, bool_Hopf)
%
% firs_Hopf_data    info on the two Hopf bifurcations:
% 
% firs_Hopf_data = cell(2)
% each cell a structure of the form x, par_alpha, par_beta, eigenvec,
% eigenval (where the eigenpair is complex)
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

new_nodes = 2:7;
new_simplices = 1:6;
list_of_Hopf_data = cell(1);
list_of_Hopf_data{1} = first_Hopf_data;
[list_of_nodes, list_of_simplices, list_of_Hopf_data] = ...
    Hopf_update_scal_eqs(list_of_nodes,...
    list_of_simplices, new_nodes, new_simplices, list_of_Hopf_data, rhs, phi);

plot(list_of_simplices);

% add the appropriate continuation equations
for k = 1:6
    simplex_k = list_of_simplices.simplex{k};
    node_numbers = simplex_k.nodes_number;
    for j = node_numbers
        % n_near_nodes = neighboring_nodes(j, list_of_nodes, list_of_simplices);
        [list_of_nodes] = simplex_scalar_equations(j, ...
            bool_Hopf, list_of_nodes);
    end
end


% validate all new simplices
for j = 1:-6
    simplex_j = list_of_simplices.simplex{j};
    use_intlab = 1;
    [flag,Imin,Imax,Yvector,Z0vector,Z1vector,Z2vector,...
        simplex_j, list_of_nodes] = ...
        radii_polynomials_simplex(simplex_j, list_of_nodes);
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
    % storage
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
    
    
    new_nodes = list_of_new_frontal_nodes(list_of_new_frontal_nodes>0);
    [list_of_nodes, list_of_simplices, list_of_Hopf_data] = ...
        Hopf_update_scal_eqs(list_of_nodes, list_of_simplices,...
        new_nodes, index_new_simplices, list_of_Hopf_data, rhs, phi);
    
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


function [list_of_nodes, list_of_simplices,list_of_Hopf_data] = ...
    Hopf_update_scal_eqs...
    (list_of_nodes,list_of_simplices, new_nodes, ...
    index_new_simplices,list_of_Hopf_data, rhs, phi)

for i = 1:length(index_new_simplices)
    simplex_i = list_of_simplices.simplex{index_new_simplices(i)};
    new_nodes_simplex = intersect(new_nodes, simplex_i.nodes_number);
    [simplex_i,list_of_Hopf_data] = Hopf_update_scal_eqs_simplex...
        (simplex_i,list_of_Hopf_data,new_nodes_simplex, rhs, phi);
    list_of_simplices.simplex{index_new_simplices(i)} = simplex_i;
    for j = 1:length(new_nodes_simplex)
        list_of_nodes{new_nodes_simplex(j)} = ...
            simplex_i.nodes{simplex_i.nodes_number==new_nodes_simplex(j)};
    end
    new_nodes = setdiff(new_nodes, new_nodes_simplex); 
    % remove already considered nodes
end

end


function [simplex_i,list_of_Hopf_data] = Hopf_update_scal_eqs_simplex...
    (simplex_i, list_of_Hopf_data, new_nodes_simplex, rhs, phi)
if length(new_nodes_simplex)==3
    error('No initial node to learn from!')
elseif length(new_nodes_simplex)==2 && simplex_i.number>1
    warning('Is this okay?')
end

old_nodes = setdiff(simplex_i.nodes_number, new_nodes_simplex);
old_node = old_nodes(1);

% we just need one info about where the two Hopf bifs are
Hopf_data = list_of_Hopf_data{old_node};

list_of_Hopf_data{end+1} = cell(2,1);

for i = 1: length(new_nodes_simplex)
    node_i = simplex_i.nodes{i};
    parameter = node_i.solution.scalar(3);
    
    % compute the two Hopf bifurcations and store the data
    % Hopf_data stores twice x, par_alpha, par_beta, eigenvec, eigenval
    x1 = Hopf_data{1}.x;
    par_alpha1 = Hopf_data{1}.par_alpha;
    eigenvec1 = Hopf_data{1}.eigenvec;
    eigenval1 = Hopf_data{1}.eigenval;
    [x1, par_alpha1,  eigenvec1, eigenval1] = Hopf_numerical(x1, par_alpha1, parameter, eigenvec1, eigenval1, rhs, phi);
    
    x1_Hopf_data.x = x1;
    x1_Hopf_data.par_alpha = par_alpha1;
    x1_Hopf_data.eigenvec = eigenvec1;
    x1_Hopf_data.eigenval = eigenval1;
    x1_Hopf_data.par_beta = parameter;
    list_of_Hopf_data{end}{1} = x1_Hopf_data;
    
    x2 = Hopf_data{2}.x;
    par_alpha2 = Hopf_data{2}.par_alpha;
    eigenvec2 = Hopf_data{2}.eigenvec;
    eigenval2 = Hopf_data{2}.eigenval;
    [x2, par_alpha2,  eigenvec2, eigenval2] = Hopf_numerical(x2, par_alpha2, parameter, eigenvec2, eigenval2,rhs, phi);
    
    x2_Hopf_data.x = x2;
    x2_Hopf_data.par_alpha = par_alpha2;
    x2_Hopf_data.eigenvec = eigenvec2;
    x2_Hopf_data.eigenval = eigenval2;
    x2_Hopf_data.par_beta = parameter;
    list_of_Hopf_data{end}{2} = x2_Hopf_data;
    
    % here we need to use the new info to update the node
    % first the problem
    node_scalar_eqs = big_Hopf_scalar_eqs(node_i.solution, list_of_Hopf_data{end});
    node_i.problem.scalar_equations = node_scalar_eqs;
    problem = node_i.problem;
    problem = continuation_equation_simplex(problem, node_i.solution);
    % then the solution
    better_sol = Newton_2(node_i.solution, problem);
    node_i.solution = better_sol;
    
    simplex_i.nodes{i} = node_i;
end

end