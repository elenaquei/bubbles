classdef simplex
    properties %(SetAccess=private)
        number
        nodes_number % list of 3 nodes
        nodes
        verified
        frontal
        verification_coeff
    end
    methods
        function simplex_x = simplex(nodes_number, number, verified, frontal, list_nodes)
            if nargin <3 || isempty(verified)
                verified = 0;
            end
            if nargin<4 || isempty(frontal)
                frontal = 1;
            end
            simplex_x.number = number;
            simplex_x.verified = verified;
            simplex_x.frontal = frontal;
            simplex_x.nodes_number = nodes_number;
            simplex_x.verification_coeff = 0; % for later, "how easy it was to verify"
            if nargin > 4
                if length(list_nodes)>3
                    list_of_nodes = list_nodes;
                    nodes = cell(3,1);
                    for i = 1:3
                        nodes{i} = list_of_nodes{nodes_number(i)};
                    end
                else 
                    nodes = list_nodes;
                end
                simplex_x.nodes = nodes;
            end
            
            
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
        
        
        function plot_simplex(simplex, color)
            % function plot_simplex(simplex, color)
            if nargin<3
                color = 'r';
            end
            x_coord = zeros(1,3);
            y_coord = zeros(1,3);
            z_coord = zeros(1,3);
            label = cell(3,1);
            for i = 1:3
                x_coord(i) = simplex.nodes{i}.solution.scalar(2);
                y_coord(i) = simplex.nodes{i}.solution.scalar(3);
                z_coord(i) = simplex.nodes{i}.solution.scalar(4);
                label{i} = 'x'+string(simplex.nodes{i}.number);
            end
            fill3(x_coord,y_coord,z_coord,color)
            alpha 0.5
            text(x_coord,y_coord,z_coord,label,'VerticalAlignment','bottom','HorizontalAlignment','right')
        end
    end
end