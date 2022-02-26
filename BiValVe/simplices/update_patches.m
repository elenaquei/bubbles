function list_of_nodes = update_patches(list_of_nodes,list_of_simplices, indices)

if length(list_of_nodes)==1
    list_of_nodes_cell = cell(1,1);
    list_of_nodes_cell{1} = list_of_nodes;
    list_of_nodes = list_of_nodes_cell;
end

if nargin <3 || isempty(indices)
    indices = 1:length(list_of_nodes);
end

for span = 1:length(indices)
    j = indices(span);
    list_of_nodes{j}.patch = [];
    for k = 1: length(list_of_simplices)
        list_of_nodes{j} = add_simplex_to_patch(list_of_nodes{j},...
            list_of_simplices.simplex{k});
    end
end
if length(list_of_nodes)==1
    list_of_nodes = list_of_nodes{1};
end
end


function node = add_simplex_to_patch(node, simplex)
if any(simplex.nodes_number == node.number)
    node.patch(end+1)= simplex.number;
end
end