function [list_of_simplices, list_of_nodes] = merge_data(save_file1, save_file2)

try 
    load save_file1 temp_list_of_nodes
    list_of_nodes1 = temp_list_of_nodes;
    clear temp_list_of_nodes
catch
    try
        load save_file1 list_of_nodes
        list_of_nodes1 = list_of_nodes;
        clear temp_list_of_nodes
    catch
        error('No node list in first file')
    end
end

load save_file1 list_of_simplices
list_of_simplices1 = list_of_simplices;
clear list_of_simplices


try 
    load save_file2 temp_list_of_nodes
    list_of_nodes2 = temp_list_of_nodes;
    clear temp_list_of_nodes
catch
    try
        load save_file2 list_of_nodes
        list_of_nodes2 = list_of_nodes;
        clear temp_list_of_nodes
    catch
        error('No node list in second file')
    end
end


load save_file2 list_of_simplices
list_of_simplices2 = list_of_simplices;
clear list_of_simplices


% merge the list of nodes
max_length_node = max( length(list_of_nodes2), length(list_of_nodes1));
list_of_nodes = cell(max_length_node,1);
for merge_i = 1: length(temp_list_of_nodes)
    if ~isempty(list_of_nodes2{merge_i})
        list_of_nodes{merge_i} = list_of_nodes2{merge_i};
    elseif ~isempty(list_of_nodes1{merge_i})
        list_of_nodes{merge_i} = list_of_nodes1{merge_i};
    end
end

% merge the list of simplices = if I was sure the ORDER was correct, I
% could use that... 
% let's try this! 
list_of_simplices = all_simplices();
count1 = 0;
count2 = 0;
for merge_i = 1: length(temp_list_of_nodes)
    if list_of_simplices1.order(count1) == merge_i
        list_of_simplices.append(list_of_simplices1.simplex{count1});
        count1 = count1 + 1;
    elseif list_of_simplices2.order(count1) == merge_i
        list_of_simplices.append(list_of_simplices2.simplex{count1});
        count2 = count2 + 1;
    else
        error('A simplex of order %i is missing',merge_i)
    end
end


end