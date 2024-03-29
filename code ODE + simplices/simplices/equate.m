function [i_loc,list_of_nodes, list_of_simplices,list_of_frontal_nodes, changed_simplices] = ...
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

new_node = node_i;
new_node.solution = coords; % this way, it maintains the patch too!
new_node.patch = [];
list_of_nodes{i_loc} = new_node;

for bigger_than_j = j_loc:length(list_of_nodes)-1
    list_of_nodes{bigger_than_j}.number = list_of_nodes{bigger_than_j}.number-1;
    list_of_nodes{bigger_than_j} = list_of_nodes{bigger_than_j+1};
end
list_of_nodes(end) = [];

changed_patches =[];
complete_patch = union(node_i.patch, node_j.patch);
for i = 1:length(complete_patch)
    index = complete_patch(i);
    simplex_i = list_of_simplices.simplex{index};
    simplex_i.nodes_number(simplex_i.nodes_number==j_loc) = i_loc;
    simplex_i.verified = 0;
    list_of_simplices.simplex{index} = simplex_i;
    changed_patches = union(changed_patches, simplex_i.nodes_number);
end

changed_simplices = complete_patch;
changed_patches(changed_patches>j_loc) = changed_patches(changed_patches>j_loc) -1;
% all simplices in the two patches have been changed

% remove all mentions of node j_loc
for index = 1:length(list_of_simplices)
    simplex_i = list_of_simplices.simplex{index};
    nodes_number = simplex_i.nodes_number;
    nodes_number(nodes_number>j_loc) = nodes_number(nodes_number>j_loc) -1;
    simplex_i.nodes_number = nodes_number;
    list_of_simplices.simplex{index}=simplex_i;
end

list_of_nodes = update_patches(list_of_nodes,list_of_simplices,changed_patches);

if nargin>4
    list_of_frontal_nodes = setdiff(list_of_frontal_nodes,j_loc);
else
    list_of_frontal_nodes = - i_loc;
end

end