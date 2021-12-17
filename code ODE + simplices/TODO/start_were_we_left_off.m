function new_name_file = start_were_we_left_off(name_file, new_niter, new_name_file)

if nargin <3 || isempty(new_name_file)
    new_name_file = ['new_', name_file];
end

load(name_file)

if ~exist('plotting_instructions', 'var') || isempty(plotting_instructions) 
    plotting_instructions = 10;
end

if bool_validated
    new_name_file = continuation_simplex_core(new_name_file, list_of_simplices,list_of_nodes,...
    new_niter, step_size, bool_Hopf, bool_validated, list_of_frontal_nodes, ...
    plotting_instructions, Interval,Z0_iter, Z1_iter, Z2_iter, Y_iter);
else
    new_name_file = continuation_simplex_core(new_name_file, list_of_simplices,list_of_nodes,...
    new_niter, step_size, bool_Hopf, bool_validated, list_of_frontal_nodes, plotting_instructions);
end