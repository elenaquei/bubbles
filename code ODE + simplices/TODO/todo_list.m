TODO = cell(2,1);
% TODO{1} = 'empty TODO folder by renaming / cleaning';
% TODO{2} = 'clean the "mains" in simplex ODE';
% TODO{3} = 'figure out best tail norms';
% TODO{3} = 'optimise caching of A and DH';
TODO{1} = 'store only the necessary files from Kevins_code';

TODO_list_length = length(TODO);
fprintf('\n what is left to do: \n')
for i = 1:TODO_list_length
    fprintf('   - %i: %s \n', i, TODO{i});
end
