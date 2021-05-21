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
    for k = 1:3
        if simplex_i.nodes{k}.number == j_loc || simplex_i.nodes{k}.number == i_loc
            simplex_i.nodes{k}=new_node;
        end
    end
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
    for k = 1:3
        if simplex_i.nodes{k}.number > j_loc
            simplex_i.nodes{k}.number = simplex_i.nodes{k}.number-1;
        end
    end
    simplex_i.nodes_number = nodes_number;
    list_of_simplices.simplex{index}=simplex_i;
end

% removal_list = [];
% changed_simplices = [];
% for index = 1:length(list_of_simplices)
%     % replace the node j by the node i whenever necessary
%     nodes_in_simplex = list_of_simplices.simplex{index}.nodes_number;
%     index_node = nodes_in_simplex == j_loc;
%     if ~isempty(index_node)
%         nodes_in_simplex(index_node) = i_loc;
%         
%         % target simplices with twice the same node
%         if length(find(nodes_in_simplex==i_loc))>1
%             removal_list(end+1) = index;
%         end
%         
%         % % replace higher valued nodes by their new numbering
%         % index_node = nodes_in_simplex >= j;
%         % nodes_in_simplex(index_node) = nodes_in_simplex(index_node)-1;
%         
%         list_of_simplices.simplex{index}.nodes = nodes_in_simplex;
%         changed_simplices(end+1)= index;
%     end
% end
% 
% for i = 1:length(removal_list)
%     index = removal_list(i);
%     for index_j = index:length(list_of_simplices)-1
%         list_of_simplices.simpleces{index_j} = list_of_simplices.simpleces{index_j+1};
%         list_of_simplices.simpleces{index_j}.number = list_of_simplices.simpleces{index_j}.number-1;
%         changed_simplices(changed_simplices == index_j) = index_j - 1;
%     end
%     list_of_simplices(end) = [];
% end

list_of_nodes = update_patches(list_of_nodes,list_of_simplices,changed_patches);

if nargin>4
    list_of_frontal_nodes = setdiff(list_of_frontal_nodes,j_loc);
else
    list_of_frontal_nodes = - i_loc;
end

end