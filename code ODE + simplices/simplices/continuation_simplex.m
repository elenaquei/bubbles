function [save_file] = continuation_simplex(x0, F,...
    n_iter, step_size, save_file, bool_Hopf, bool_validated, plotting_instructions)
% function list_of_simplices = continuation_simplex(x0, F, n_iter, h, save_file, bool_Hopf)
%
% plotting_instructions = n > 0 (DEFAULT  = 10) plots the new simplices every
% n iterations
% plotting_instructions = - n < 0 plots ALL simplices every
% n iterations

global use_intlab
global talkative
use_intlab = 0;
if nargin< 7 || isempty(bool_validated)
    bool_validated = 1;
end
if nargin < 8 || isempty(plotting_instructions)
    plotting_instructions = 10;
end

Interval = zeros(2,6);
Z0_iter = zeros(x0.size_scalar+x0.size_vector,6);
Z1_iter = zeros(x0.size_scalar+x0.size_vector,6);
Z2_iter = zeros(x0.size_scalar+x0.size_vector,6);
Y_iter = zeros(x0.size_scalar+x0.size_vector,6);

%% set up the first 6 simplices
first_node = node(1,x0,F);
[list_of_simplices,list_of_nodes,list_of_frontal_nodes] = ...
    simplices_set_up(first_node, F, step_size);


% validate all new simplices
if bool_validated
    for j = 1:6
        simplex_jj = list_of_simplices.simplex{j};
        use_intlab = 1;
        
        [flag,Imin,Imax,Yvector,Z0vector,Z1vector,Z2vector,...
            simplex_jj, list_of_nodes] = ...
            radii_polynomials_simplex(simplex_jj, list_of_nodes, bool_Hopf);
        use_intlab = 0;
        if flag < 1
            error('Validation failed')
        elseif talkative>0
            fprintf('The validation of the %i-th simplex succeeded\n', j)
            if talkative>1
                %figure
                plot(list_of_simplices, list_of_nodes, simplex_jj.number,'b');
                hold on
            end
        end
        simplex_jj.verified = 1;
        list_of_simplices.simplex{j} = simplex_jj;
        %storage
        %     Interval(:,j)  = [Imin,Imax]';
        %     Z0_iter(:,j)   = vert(Z0vector);
        %     Z1_iter(:,j)   = vert(Z1vector);
        Interval(:,j)  = [Imin,Imax]';
        Z0_iter(:,j)   = vert(Z0vector);
        Z1_iter(:,j)   = vert(Z1vector);
        Z2_iter(:,j)   = vert(Z2vector);
        Y_iter(:,j)    = vert(Yvector);
    end
end

plot(list_of_simplices, list_of_nodes);

if bool_validated
save_file = continuation_simplex_core(save_file, list_of_simplices,list_of_nodes,...
    n_iter, step_size, bool_Hopf, bool_validated, list_of_frontal_nodes, ...
    plotting_instructions, Interval,Z0_iter, Z1_iter, Z2_iter, Y_iter);
else
    save_file = continuation_simplex_core(save_file, list_of_simplices,list_of_nodes,...
    n_iter, step_size, bool_Hopf, bool_validated, list_of_frontal_nodes,...
    plotting_instructions);
    
end



