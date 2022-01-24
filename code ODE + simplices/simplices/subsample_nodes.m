function [partial_list_of_nodes, renumbered_partial_simplices] = subsample_nodes(list_of_nodes, partial_list_of_simplices)
% function partial_list_of_nodes = subsample_nodes(list_of_nodes, partial_list_of_simplices)
%
% given a list of nodes and a subsample of its associated list of
% simplices, select only the needed nodes
% the list of simplices also return with renumbered nodes, since nodes are
% defined 
partial_list_of_nodes = cell(0,0);
renumbered_partial_simplices = partial_list_of_simplices;
nodes_numbers = [];
for i = 1:length(partial_list_of_simplices)
    new_nodes = partial_list_of_simplices.simplex{i}.nodes_number;
    nodes_to_add = setdiff(new_nodes,nodes_numbers);
    nodes_numbers = [nodes_numbers, nodes_to_add];
    for j = 1:length(nodes_to_add)
        partial_list_of_nodes{end+1} = list_of_nodes{nodes_to_add(j)};
        partial_list_of_nodes{end}.number = length(partial_list_of_nodes); % renumbering of nodes!
    end
    
    % renumbering of node numbering in simplex
    new_numbering = 0*new_nodes;
    for k = 1:length(new_nodes)
        new_numbering(k) = find(nodes_numbers==new_nodes(k));
    end
    partial_list_of_simplices.simplex{i}.nodes_number = new_numbering;
end

    