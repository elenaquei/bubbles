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
        
        function [partial_list_of_simplices, partial_list_of_nodes] = ...
                subsample(list_of_simplices, indices, list_of_nodes)
            % function [partial_list_of_simplices, partial_list_of_nodes] = ...
            %    subsample(list_of_simplices, indices, list_of_nodes)
            % 
            % with two inputs, returns a list of simplices only including
            % the simplices with index in indices. With three inputs, it
            % also returns the minimum list of nodes need for the smaller
            % list of simplices and renumbers the nodes
            
            partial_list_of_simplices = all_simplices();
            partial_list_of_simplices.simplex = cell(length(indices),1);
            for j = 1:length(indices)
                partial_list_of_simplices.order(j) = list_of_simplices.order(indices(j));
                partial_list_of_simplices.simplex{j} = list_of_simplices.simplex{indices(j)} ;
            end
            if nargin>2
                [partial_list_of_nodes, renumbered_partial_simplices] = ...
                    subsample_nodes(list_of_nodes, partial_list_of_simplices);
                partial_list_of_simplices = renumbered_partial_simplices;
            end
        end
        
        function list_of_simplices = append(list_of_simplices, simplex, varargin)
            % append as many simplices as you like, one by one to the list
            list_of_simplices.order(end+1) = simplex.number;
            list_of_simplices.simplex{end+1} = simplex;
            if ~isempty(varargin)
                list_of_simplices = append(list_of_simplices, varargin{:});
            end
        end
        
        function plot(list_simplex, list_of_nodes, index_simplices, varargin)
            global talkative
            
            if nargin <3 || isempty(index_simplices)
                index_simplices = 1:length(list_simplex);
            end
            
            fail = 0;
            for index = 1: length(index_simplices)
                i = index_simplices(index);
                try
                plot_simplex(list_simplex.simplex{i}, list_of_nodes, varargin{:});
                hold on
                catch
                    if talkative && fail < 1
                        disp('The plotting of one or more simplices failed - likely because the nodes are not in memory')
                        fail = 2;
                    end
                end
                
            end
            alpha 0.5
            set(gca,'FontSize',18)
            hold off
        end
        
        
    end
end