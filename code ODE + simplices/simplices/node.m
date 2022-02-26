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
            U = make_it_real_cols(U, node.solution);
            
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
            
            edge1 = U.' * make_it_real_cols(edge1_in_Rn, node.solution);
            edge2 = U.' * make_it_real_cols(edge2_in_Rn, node.solution);
            internal_dir = U.' * make_it_real_cols(internal_coord, node.solution);
            
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



function x_real_vec = make_it_real(x_vec, x0)
% function x_real_vec = make_it_real(vec, x0)
% takes a complex vector, rotates it and turn it into a sin-cos series

% ritation to symmetry
%x_vec = vec2Xi_vec(vec,x0);
[~,index] = max(abs(real(x_vec(1:x0.size_scalar))));
angle = atan( imag(x_vec(index))/real(x_vec(index)));
x_vec = exp( - 1i * angle) * x_vec; 
% bringing x_Xi_vec to be symmetric (by multiplication with the appropriate complex rotation)
x_Xi_vec = symmetrise(vec2Xi_vec(x_vec,x0));

x_real_vec = 0*x_vec;
x_real_vec(1:x0.size_scalar) = x_Xi_vec.scalar;
for j = 1: x0.size_vector
    index_modes = x0.size_scalar + (j-1)*(2*x0.nodes+1) + (1:2*x0.nodes+1);
    index_positive_modes = index_modes(x0.nodes+1:end);
    index_negative_modes = index_modes(1:x0.nodes);
    
    real_part = real(x_Xi_vec.vector(j,x0.nodes+1:end));
    imaginary_part = -imag(x_Xi_vec.vector(j,1:x0.nodes));
    
    x_real_vec(index_negative_modes) = imaginary_part;
    x_real_vec(index_positive_modes) = real_part;
end

end


function real_cols = make_it_real_cols(cols, x0)
real_cols = 0*cols;
for i = 1: size(cols,2)
    real_cols(:,i) = make_it_real(cols(:,i),x0);
end

end

function x_complex_vec = make_it_complex(vec,x0)
% function x_real_vec = make_it_complex(vec, x0)
% takes a real vector that stores a sin-cos series and returns a "proper"
% Fourier series vector

x_complex_vec = 0*vec;
x_complex_vec(1:x0.size_scalar) = vec(1:x0.size_scalar);
for j = 1: x0.size_vector
    index_modes = x0.size_scalar + (j-1)*(2*x0.nodes+1) + (1:(2*x0.nodes+1));
    index_positive_modes = index_modes(x0.nodes+2:end);
    index_negative_modes = index_modes(x0.nodes:-1:1);
    zeroth_mode = index_modes(x0.nodes+1);
    
    real_part = vec(index_positive_modes);
    imaginary_part = vec(index_negative_modes);
    
    x_complex_vec(zeroth_mode) = vec(zeroth_mode);
    
    x_complex_vec(index_negative_modes) = real_part - 1i*imaginary_part;
    x_complex_vec(index_positive_modes) = real_part +  1i*imaginary_part;
end


% testing
x_Xi = vec2Xi_vec(x_complex_vec, x0);
if any(norm(x_Xi - symmetrise(x_Xi)) > 10^-7)
    error('this is not good')
end


end

