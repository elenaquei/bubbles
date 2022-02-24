function stored_index = store_center(save_file, save_file_temp, list_of_frontal_nodes)

load(save_file);

list_of_non_frontal_nodes = setdiff(1:length(list_of_nodes),list_of_frontal_nodes);
loc_list_of_nodes = cell(1, length(list_of_non_frontal_nodes));
copy_index = list_of_non_frontal_nodes;
%deleted_index = [];
stored_index = [];
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
            stored_index(end+1) = index;
        else
            copy_index = setdiff(copy_index, index);
        end
    else
        copy_index = setdiff(copy_index, index);
    end
end

try
    load(save_file_temp)
catch ME
    if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
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
save(save_file_temp,'temp_list_of_nodes')


end