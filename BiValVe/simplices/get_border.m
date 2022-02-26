function [border_list_of_nodes, border_list_of_simplices] = get_border(list_of_nodes, list_of_simplices)
% function [border_list_of_nodes, border_list_of_simplices] = get_border(list_of_nodes, list_of_simplices)
%
% given a simplex manifold, only the frontal simplices and the needed nodes
% are returned

non_frontal_simplices = [];

for i = 1: length(list_of_simplices)
    if ~list_of_simplices.simplex{i}.frontal
        non_frontal_simplices(end+1) = i;
    end
end

[border_list_of_nodes, border_list_of_simplices] = subsample(list_of_simplices,...
    non_frontal_simplices, list_of_nodes);

end