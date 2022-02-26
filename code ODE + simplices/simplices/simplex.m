classdef simplex
    properties %(SetAccess=private)
        number
        nodes_number % list of 3 nodes
        % nodes
        verified
        frontal
        verification_coeff
        validation_interval
    end
    methods
        function simplex_x = simplex(nodes_number, number, verified, frontal, validation_interval)
            if nargin <3 || isempty(verified)
                verified = 0;
            end
            if nargin<4 || isempty(frontal)
                frontal = 1;
            end
            if nargin < 5 || verified == 0
                validation_interval = [NaN, NaN];
            end
            simplex_x.number = number;
            simplex_x.verified = verified;
            simplex_x.frontal = frontal;
            simplex_x.nodes_number = nodes_number;
            simplex_x.verification_coeff = 0; % for later, "how easy it was to verify"
            simplex_x.validation_interval = validation_interval;
        end
        
        function loc_node = node(simplex, list_of_nodes, i)
            node_number = simplex.nodes_number(i);
            loc_node = list_of_nodes{node_number};
        end
        
        
        function [simplex_x,Imin,Imax,Yvector,Z0vector,Z1vector,Z2vector,Ys]...
                = verify(simplex_x, varargin)
            % function [simplex_x,Imin,Imax,Yvector,Z0vector,Z1vector,Z2vector,Ys]...
            %    = verify(simplex_x, varargin)
            % additional arguments in varargin
            % DH0,DH1,DH2,A0,A1,A2
            
            xBar0 = simplex_x.nodes{1}.solution;
            xBar1 = simplex_x.nodes{2}.solution;
            xBar2 = simplex_x.nodes{3}.solution;
            alpha0 = simplex_x.nodes{1}.problem;
            alpha1 = simplex_x.nodes{2}.problem;
            alpha2 = simplex_x.nodes{3}.problem;
            previous_iter0 = simplex_x.nodes{1}.previous_validation;
            previous_iter1 = simplex_x.nodes{2}.previous_validation;
            
            [flag,Imin,Imax,new_iter,Yvector,Z0vector,Z1vector,Z2vector,new_step,Ys] = ...
                radii_polynomials_simplex(xBar0,xBar1,xBar2,...
                alpha0,alpha1,alpha2,previous_iter0,previous_iter1,varargin{:});
            simplex_x.verified = flag;
            simplex_x.verification_coeff = new_step;
            simplex_x.node{3}.previous_validation = new_iter;
        end
        
        
        function plot_simplex(simplex, list_of_nodes, color, variables, edge_alpha, plot_sums)
            % function plot_simplex(simplex, color)
            if nargin<3
                if simplex.verified
                    color = 'b';
                elseif simplex.frontal
                    color = 'y';
                else
                    color = 'r';
                end
            end
            if nargin<4
                variables = [2,3,4];
            end
            if nargin<5
                edge_alpha = 1;
            end
            if nargin<6
                plot_sums = [NaN;NaN;NaN];
            end
            x_coord = zeros(1,3);
            y_coord = zeros(1,3);
            z_coord = zeros(1,3);
            label = cell(3,1);
            for i = 1:3
                if isnan(plot_sums(1))
                    x_coord(i) = list_of_nodes{simplex.nodes_number(i)}.solution.scalar(variables(1));
                else
                    x_coord(i) = real(sum(list_of_nodes{simplex.nodes_number(i)}.solution.vector(plot_sums(1),:)));
                end
                if isnan(plot_sums(2))
                    y_coord(i) = list_of_nodes{simplex.nodes_number(i)}.solution.scalar(variables(2));
                else
                    y_coord(i) = real(sum(list_of_nodes{simplex.nodes_number(i)}.solution.vector(plot_sums(2),:)));
                end
                if isnan(plot_sums(3))
                    z_coord(i) = list_of_nodes{simplex.nodes_number(i)}.solution.scalar(variables(3));
                else
                    z_coord(i) = real(sum(list_of_nodes{simplex.nodes_number(i)}.solution.vector(plot_sums(3),:)));
                end
                %x_coord(i) = list_of_nodes{simplex.nodes_number(i)}.solution.scalar(6);
                %y_coord(i) = list_of_nodes{simplex.nodes_number(i)}.solution.scalar(5);
                %z_coord(i) = list_of_nodes{simplex.nodes_number(i)}.solution.scalar(3);
                label{i} = 'x'+string(simplex.nodes_number(i));
            end
            label_simplex = 'S' +string(simplex.number);
            center_x = mean(x_coord);
            center_y = mean(y_coord);
            center_z = mean(z_coord);
            fill3(x_coord,y_coord,z_coord,color,'EdgeAlpha',edge_alpha)
            
            %text(x_coord,y_coord,z_coord,label,'VerticalAlignment','bottom','HorizontalAlignment','right')
            %text(center_x,center_y,center_z,label_simplex,'VerticalAlignment','bottom','HorizontalAlignment','right')
        end
    end
end
