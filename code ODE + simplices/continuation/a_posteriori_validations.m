function [list_of_simplices, index_failed_validation, Interval, Z0_iter, ...
    Z1_iter, Z2_iter, Y_iter] = a_posteriori_validations(list_of_simplices,...
    list_of_nodes, index_to_validate_simplices, bool_Hopf, max_refinements, save_file)
global talkative

if talkative>1
    fprintf('\nStarting validation of simplices, time %s\n\n',datestr(now,13));
end
if nargin == 1
    load(list_of_simplices);
end

if nargin < 4 || isempty(bool_Hopf)
    bool_Hopf =0;
end
if nargin < 3 || isempty(index_to_validate_simplices)
    index_to_validate_simplices = 1:length(list_of_simplices);
end
if nargin < 5 || isempty(max_refinements)
    max_refinements = 5;
end
if nargin < 6
    save_file = 'validation_results';
end

ref = 0;
index_failed_validation = [];
Interval = [];
Z0_iter = [];
Z1_iter = [];
Z2_iter = [];
Y_iter= [];
temp_save_file = append('temp_', save_file);

while ~isempty(index_to_validate_simplices)
    
    [list_of_simplices, index_failed_validation, Interval, Z0_iter, ...
        Z1_iter, Z2_iter, Y_iter] = iteration_validation(list_of_simplices,...
        list_of_nodes, index_to_validate_simplices, bool_Hopf, temp_save_file);
    ref = ref + 1;
    save(save_file,'list_of_simplices','list_of_nodes','Interval',...
        'Z0_iter','Z1_iter','Z2_iter','Y_iter','bool_Hopf');
    if ref > max_refinements
        break
    end
    if ~isempty(index_failed_validation)
        [list_of_simplices, list_of_nodes, index_to_validate_simplices] = ...
            refine_simplices(list_of_simplices,...
            list_of_nodes, index_failed_validation);
    else
        break
    end
end
end

function [list_of_simplices, index_non_validated, Interval, Z0_iter, ...
    Z1_iter, Z2_iter, Y_iter] = iteration_validation(list_of_simplices,...
    list_of_nodes, index_validation_simplices, bool_Hopf, save_file)
global use_intlab
global talkative

n_iter = length(index_validation_simplices);

Interval = zeros(2,n_iter);
x0 = list_of_nodes{1}.solution;
Z0_iter = zeros(x0.size_scalar+x0.size_vector,n_iter);
Z1_iter = zeros(x0.size_scalar+x0.size_vector,n_iter);
Z2_iter = zeros(x0.size_scalar+x0.size_vector,n_iter);
Y_iter = zeros(x0.size_scalar+x0.size_vector,n_iter);
if nargin > 7
    Interval = [ Interval_old, Interval];
    Z0_iter = [Z0_iter_old, Z0_iter];
    Z1_iter = [Z1_iter_old, Z1_iter];
    Z2_iter = [Z2_iter_old, Z2_iter];
    Y_iter = [Y_iter_old, Y_iter];
end
index_non_validated = [];

for index_j = 1:length(index_validation_simplices)
    j = index_validation_simplices(index_j);
    simplex_jj = list_of_simplices.simplex{j};
    if simplex_jj.verified
        continue
    end
    try
        use_intlab = 1;
        [flag,Imin,Imax,Yvector,Z0vector,Z1vector,Z2vector,...
            simplex_jj, list_of_nodes] = ...
            radii_polynomials_simplex(simplex_jj, list_of_nodes,bool_Hopf);
        use_intlab = 0;
        
        if flag < 1
            index_non_validated(end+1) = j;
            if talkative > 1
                fprintf('The validation of the %i-th simplex failed\n', j)
            end
        elseif talkative>0
            fprintf('The validation of the %i-th simplex succeeded\n', j)
        end
        if flag
            simplex_jj.verified = 1;
            simplex_jj.validation_interval = [Imin, Imax];
            list_of_simplices.simplex{j} = simplex_jj;
            %storage
            Interval(:,j)  = [Imin,Imax]';
            Z0_iter(:,j)   = vert(Z0vector);
            Z1_iter(:,j)   = vert(Z1vector);
            Z2_iter(:,j)   = vert(Z2vector);
            Y_iter(:,j)    = vert(Yvector);
        end
    catch ME
        if (strcmp(ME.identifier, 'VAL:NegDelta'))
            use_intlab = 0;
            index_non_validated(end+1) = j;
        else
            rethrow(ME)
        end
    end
    
    if any(mod(index_j,5)==0)
        save(save_file,'list_of_simplices','list_of_nodes','Interval',...
            'Z0_iter','Z1_iter','Z2_iter','Y_iter','bool_Hopf');
    end
end
end


function [list_of_simplices, list_of_nodes, index_new_simplices] = ...
    refine_simplices(list_of_simplices,...
    list_of_nodes, index_to_refine)

index_new_simplices = [];

for j = index_to_refine
    old_simplex = list_of_simplices.simplex{j};
    frontal = old_simplex.frontal;
    old_node1 = list_of_nodes{old_simplex.nodes_number(1)};
    old_node2 = list_of_nodes{old_simplex.nodes_number(2)};
    old_node3 = list_of_nodes{old_simplex.nodes_number(3)};
    
    node1_number = old_node1.number;
    node2_number = old_node2.number;
    node3_number = old_node3.number;
    node12_number = length(list_of_nodes)+1;
    node23_number = length(list_of_nodes)+2;
    node13_number = length(list_of_nodes)+3;
    
    new_node12 = mean_node(old_node1, old_node2, node12_number);
    new_node23 = mean_node(old_node3, old_node2, node23_number);
    new_node13 = mean_node(old_node1, old_node3, node13_number);
    
    list_of_nodes{node12_number} = new_node12;
    list_of_nodes{node23_number} = new_node23;
    list_of_nodes{node13_number} = new_node13;
    
    number_new_simplex1 = j;
    number_new_simplex2 = length(list_of_simplices)+1;
    number_new_simplex3 = length(list_of_simplices)+2;
    number_new_simplex4 = length(list_of_simplices)+3;
    numbers_new_simplices = [number_new_simplex1, number_new_simplex2, number_new_simplex3, number_new_simplex4];
    
    new_simplex1 = simplex([node1_number, node12_number, node13_number], number_new_simplex1,0,frontal);
    new_simplex2 = simplex([node2_number, node12_number, node23_number], number_new_simplex2,0,frontal);
    new_simplex3 = simplex([node3_number, node13_number, node23_number], number_new_simplex3,0,frontal);
    new_simplex4 = simplex([node23_number, node12_number, node13_number], number_new_simplex4,0,frontal);
    
    list_of_simplices.simplex{j} = new_simplex1;
    list_of_simplices = append(list_of_simplices, new_simplex2, new_simplex3, new_simplex4);
    index_new_simplices(end+1:end+4) = numbers_new_simplices;
end

end


function new_node = mean_node(old_node1, old_node2, node_number)
solution1 = old_node1.solution;
solution2 = old_node2.solution;
mean_solution = 0.5*solution1 + 0.5*solution2;

problem1 = old_node1.problem;
problem2 = old_node2.problem;

mean_problem = problem1;
for i = 1:3
    mean_problem.scalar_equations.linear_coef{i} = ...
        0.5*(problem1.scalar_equations.linear_coef{i}+ ...
        problem2.scalar_equations.linear_coef{i});
end
new_node = node(node_number, mean_solution, mean_problem);
end
