classdef node
    properties
        number
        solution
        problem
        previous_validation
        patch
    end
    methods
        function node_x = node(number, solution, problem, patch, previous_val)
            % solution is a Xi_vector
            % problem is a full_problem
            % previous_val is a struct with Y, Z0, Z1 vectors
            if nargin < 4 || isempty(patch)
                patch = [];
            end
            if nargin < 5 || isempty(previous_val)
                previous_val = [];
            end
            node_x.number = number;
            node_x.solution = solution;
            node_x.problem = problem;
            node_x.patch = patch;
            node_x.previous_validation = previous_val;
        end
        
        function coord = coordinate(node,i)
            coord = node.solution.scalar(i+1);
        end
        
        function plot_patch(node,list_simplex)
            for i = 1: length(node.patch)
                plot_simplex(list_simplex.simplex{node.patch(i)},'b');
                hold on
            end
            hold off
        end
        
        function [out_nodes,all_nodes] = open_edges(node, list_of_simplices)
            % function [out_nodes,all_nodes] = open_edges(node, list_of_simplices)
            out_nodes = [];
            all_nodes = [];
            for i = node.patch
                for index_j = 1:length(list_of_simplices.simplex{i}.nodes_number)
                    j = list_of_simplices.simplex{i}.nodes_number(index_j);
                    if j == node.number
                        continue
                    end
                    if any(out_nodes == j)
                        out_nodes(out_nodes==j) = []; % if it's in, cancel it
                    else
                        out_nodes(end+1) = j;
                    end
                    all_nodes(end+1) = j;
                end
            end
            all_nodes = unique(all_nodes);
        end
        
        
        
        function [gap] = gap_angle(node, DF_x, list_of_simplices, list_of_nodes)
            
            U = null(DF_x);
            
            [out_nodes,all_nodes] = open_edges(node, list_of_simplices);
            
            if length(out_nodes)~=2
                error('there are not two outnodes?!')
            end
            internal_nodes = setdiff(all_nodes,out_nodes);
            if isempty(internal_nodes)
                internal_coord = 1/2 *(list_of_nodes{out_nodes(1)}.coordinates + ...
                    list_of_nodes{out_nodes(2)}.coordinates); % average of two open edges
            else
                internal_coord = list_of_nodes{internal_nodes(1)}.solution-node.solution;
            end
            internal_coord = Xi_vec2vec(internal_coord) / max(norm(internal_coord));
            % open edges
            solution_vec = Xi_vec2vec(node.solution);
            edge1_in_Rn = Xi_vec2vec(list_of_nodes{out_nodes(1)}.solution);
            edge1_in_Rn = (edge1_in_Rn - solution_vec) / norm(edge1_in_Rn);
            edge2_in_Rn = Xi_vec2vec(list_of_nodes{out_nodes(2)}.solution);
            edge2_in_Rn = (edge2_in_Rn - solution_vec) / norm(edge2_in_Rn);
            
            if size(edge1_in_Rn,1) == 1
                edge1_in_Rn=edge1_in_Rn.';
            end
            if size(edge2_in_Rn,1) == 1
                edge2_in_Rn=edge2_in_Rn.';
            end
            if size(internal_coord,1) == 1
                internal_coord=internal_coord.';
            end
            
            edge1 = U.' * edge1_in_Rn;
            edge2 = U.' * edge2_in_Rn;
            internal_dir = U.' * internal_coord;
            
            rotation_angle = angle(edge2);
            rotation_matrix = [ cos(rotation_angle), sin(rotation_angle);
                -sin(rotation_angle), cos(rotation_angle)];
            
            rotated_edge1 = rotation_matrix * edge1;
            rotated_internal_dir = rotation_matrix * internal_dir;
            
            angle_edge1 = angle(rotated_edge1);
            angle_internal_dir = angle(rotated_internal_dir);
            
            if angle_edge1 > angle_internal_dir
                gap = -(2 * pi - angle_edge1);
            else
                gap = angle_edge1;
            end
            
            
        end
    end
end
