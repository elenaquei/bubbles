function [i_loc,list_of_nodes, list_of_simplices,list_of_frontal_nodes] = ...
    equate(list_of_nodes, i_loc,j_loc, ...
    list_of_simplices, list_of_frontal_nodes)
% replace nodes i and j with their average, and return the number of the
% new node

if i_loc>j_loc
    temp = i_loc;
    i_loc = j_loc;
    j_loc = temp;
end

node_i = list_of_nodes{i_loc};
node_j = list_of_nodes{j_loc};

coords = 1/2 * (node_i.solution + node_j.solution);

if node_i.frontal ~= node_j.frontal
    warning('Merging of a frontal node with a non-frontal node')
end
new_node = node_i;
new_node.solution = coords; % this way, it maintains the patch too!
new_node.patch = [];
list_of_nodes{i_loc} = new_node;

for bigger_than_j = j_loc:length(list_of_nodes)-1
    list_of_nodes{bigger_than_j} = list_of_nodes{bigger_than_j+1};
    
end
list_of_nodes(end) = [];

removal_list = [];
for index = 1:length(list_of_simplices)
    % replace the node j by the node i whenever necessary
    nodes_in_simplex = list_of_simplices.simplices{index}.nodes;
    index_node = nodes_in_simplex == j_loc;
    nodes_in_simplex(index_node) = i_loc;
    
    % target simplices with twice the same node
    if length(find(nodes_in_simplex==i_loc))>1
        removal_list(end+1) = index;
    end
    
    % % replace higher valued nodes by their new numbering
    % index_node = nodes_in_simplex >= j;
    % nodes_in_simplex(index_node) = nodes_in_simplex(index_node)-1;
    
    list_of_simplices.simpleces{index}.nodes = nodes_in_simplex;
end

for index = removal_list
    for index_j = index:length(list_of_simplices)-1
        list_of_simplices.simpleces{index_j} = list_of_simplices.simpleces{index_j+1};
        list_of_simplices.simpleces{index_j}.number = list_of_simplices.simpleces{index_j}.number-1;
    end
    list_of_simplices(end) = [];
end

list_of_nodes = update_patches(list_of_nodes,list_of_simplices);

if nargin>4
    list_of_frontal_nodes = setdiff(list_of_frontal_nodes,j_loc);
end

end