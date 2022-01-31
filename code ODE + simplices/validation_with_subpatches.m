function validation_with_subpatches(save_file, save_file_iter, subsections)
global talkative

load(save_file)

if nargin<3
    subsections = 200;
end
disp('Subpatches validation started')
disp('')
size_subsections = floor(length(list_of_simplices)/subsections);

if talkative
    fprintf('There are %i subpatches of size %i\n', subsections,size_subsections)
end

n_simplices = length(list_of_simplices);
last_index = 0;
for i = 1:subsections
    if i == subsections
        indices = last_index+1:n_simplices;
    else
        indices = last_index+(1:size_subsections);
    end
    if talkative
        fprintf("Validation of patch %i out of %i started \n ", i,subsections)
    end
    talkative = talkative-1;
    last_index = last_index + size_subsections;
    load(save_file)
    [partial_list_of_simplices, partial_list_of_nodes] = ...
                subsample(list_of_simplices, indices, list_of_nodes);
    clear list_of_simplices list_of_nodes
    save_file_iteri = append(save_file_iter,num2str(i));
    a_posteriori_validations(partial_list_of_simplices,...
        partial_list_of_nodes, [], bool_Hopf,[],save_file_iteri);
    talkative = talkative+1;
    if talkative 
        fprintf("Validation of patch %i finished \n", i)
    end
end

end