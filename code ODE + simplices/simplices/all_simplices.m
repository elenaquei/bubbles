classdef all_simplices 
    properties
        order
        simplex
    end
    methods
        function list_of_simplices = all_simplices()
            % create an empty list of nodes
            list_of_simplices.order =[];
            list_of_simplices.simplex = cell(0,0);
        end
        
        function len = length(list_of_simplices)
            len = length(list_of_simplices.order);
        end
        
        function list_of_simplices = append(list_of_simplices, simplex, varargin)
            % append as many simplices as you like, one by one to the list
            list_of_simplices.order(end+1) = simplex.number;
            list_of_simplices.simplex{end+1} = simplex;
            if ~isempty(varargin)
                list_of_simplices = append(list_of_simplices, varargin{:});
            end
        end
        
        function plot(list_simplex)
            for i = 1: length(list_simplex)
                plot_simplex(list_simplex.simplex{i});
                hold on
            end
            hold off
        end
        
        
    end
end