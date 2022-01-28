function save_file = continuation_simplex_core(save_file, list_of_simplices,list_of_nodes,...
    n_iter, step_size, bool_Hopf, bool_validated, list_of_frontal_nodes, ...
    plotting_instructions, Interval_old,Z0_iter_old, Z1_iter_old, Z2_iter_old, Y_iter_old)
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
use_intlab = 0;
F = list_of_nodes{end}.problem;

save_file_nodes = append(save_file,'_nodes_temp.mat');
prompt = sprintf('If it exists, the file %s will be deleted, press 0 to avoid this, but program will likely crash.\nIf no info is given within 3 seconds, the computation will continue automatically', save_file_nodes);
t = timer('ExecutionMode', 'singleShot', 'StartDelay', 3, 'TimerFcn', @pressEnter);
start(t);
bool = input(prompt);
stop(t);
delete(t);

if isempty(bool) || bool~=0
    try
        delete(save_file_nodes)
    catch
    end
end

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
            simplex_jj.validation_interval = [Imin, Imax];
            list_of_simplices.simplex{j} = simplex_jj;
            %storage
            Interval(:,j)  = [Imin,Imax]';
            Z0_iter(:,j)   = vert(Z0vector);
            Z1_iter(:,j)   = vert(Z1vector);
            Z2_iter(:,j)   = vert(Z2vector);
            Y_iter(:,j)    = vert(Yvector);
        end
    end
    
    if any(mod(i:i+length(index_new_simplices),plotting_instructions)==0)
        hold on
        if plotting_instructions >0
            plot_index = union(plot_index,index_new_simplices);
            plot(list_of_simplices, list_of_nodes,plot_index);
        else
            plot(list_of_simplices, list_of_nodes)
        end
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
    
    
    %if talkative>1
    %    figure
    %    plot(list_of_simplices,list_of_nodes);
    %end
    i = i + length(index_new_simplices);
    
    
    % memory management: if there are many non-frontal modes, store them in
    % a file and "forget" them
    if any(mod(i:i+length(index_new_simplices),250)==0)
        list_of_non_frontal_nodes = setdiff(1:length(list_of_nodes),list_of_frontal_nodes);
        loc_list_of_nodes = cell(1, length(list_of_non_frontal_nodes));
        copy_index = list_of_non_frontal_nodes;
        deleted_index = [];
        for delete_i = 1:length(list_of_non_frontal_nodes)
            index = list_of_non_frontal_nodes(delete_i);
            if ~isempty(list_of_nodes{index})
                loc_patch = list_of_nodes{index}.patch;
                
                all_neighbours = [];
                for check_delete_i = 1: length(loc_patch)
                    simplex_di = list_of_simplices.simplex{loc_patch(check_delete_i)};
                    all_neighbours = union(all_neighbours, simplex_di.nodes_number);
                end
                if ~isempty(intersect(all_neighbours, list_of_frontal_nodes))
                    bool_delete = 0;
                else
                    bool_delete = 1;
                end
                if bool_delete
                    loc_list_of_nodes{index} = list_of_nodes{index};
                    deleted_index(end+1) = index;
                else
                    copy_index = setdiff(copy_index, index);
                end
            else
                copy_index = setdiff(copy_index, index);
            end
        end
        
        try
            load(save_file_nodes)
        catch ME
            if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile') && i < 500
                temp_list_of_nodes = cell(0,0);
            else
                rethrow(ME)
            end
        end
        % merge loc_list_of_nodes and temp_list_of_nodes
        temp_list_of_nodes{length(temp_list_of_nodes)+length(copy_index)} = []; % preallocation
        for merge_i =  1:length(copy_index)
            temp_list_of_nodes{copy_index(merge_i)} = loc_list_of_nodes{copy_index(merge_i)};
        end
        save(save_file_nodes,'temp_list_of_nodes')
        clear temp_list_of_nodes
        
        for delete_i = 1:length(deleted_index)
            index = deleted_index(delete_i);
            list_of_nodes{index}=[];
        end
    end
    
end

% merge all info on nodes back again
try
    load(save_file_nodes)
catch
    temp_list_of_nodes = cell(0,0);
end
% merge loc_list_of_nodes and temp_list_of_nodes
for merge_i = 1: length(temp_list_of_nodes)
    list_of_nodes{merge_i} = temp_list_of_nodes{merge_i};
end
delete save_file_nodes

if bool_validated
save(save_file,'list_of_simplices','list_of_nodes','Interval','Z0_iter',...
    'Z1_iter','Z2_iter','Y_iter','step_size','bool_Hopf', 'bool_validated',...
    'list_of_frontal_nodes');
else
    save(save_file,'list_of_simplices','list_of_nodes',...
    'step_size','bool_Hopf', 'bool_validated',...
    'list_of_frontal_nodes');
end

if talkative>0
    fprintf('The growth is completed\n')
end

end

function priority_node = select_growing_node(list_of_simplices, list_of_nodes, list_of_frontal_nodes, F)
gap_angle_vec = 0 * list_of_frontal_nodes;
epsilon_vec = 0 * list_of_frontal_nodes;
a_vec = 0 * list_of_frontal_nodes;
for i=1:length(list_of_frontal_nodes)
    node_index = list_of_frontal_nodes(i);
    node = list_of_nodes{node_index};
    DF_x = derivative_to_matrix(derivative(F, node.solution,0));
    gap_angle_vec(i) = gap_angle(node, DF_x, list_of_simplices, list_of_nodes);
    epsilon_vec(i) = node.solution.scalar(3);
    a_vec(i) = node.solution.scalar(4);
end
if any(a_vec>=-2E-2)
    gap_angle_vec = gap_angle_vec(a_vec>=-2E-2);
    epsilon_vec = epsilon_vec(a_vec>=-2E-2);
    list_of_frontal_nodes = list_of_frontal_nodes(a_vec>=-2E-2);
end


[~, priority_node_index] = min(abs(gap_angle_vec));
min_gap = gap_angle_vec(priority_node_index);
if abs(min_gap) > pi/5
    priority_node = min(list_of_frontal_nodes);
    if -priority_node > 100
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


function pressEnter(HObj, event)
  import java.awt.*;
  import java.awt.event.*;
  rob = Robot;
  rob.keyPress(KeyEvent.VK_ENTER)
  rob.keyRelease(KeyEvent.VK_ENTER)
end