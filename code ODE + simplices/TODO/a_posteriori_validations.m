function [list_of_simplices, index_non_validated, Interval, Z0_iter, ...
    Z1_iter, Z2_iter, Y_iter] = a_posteriori_validations(list_of_simplices,...
    list_of_nodes, index_validation_simplices, bool_Hopf, max_refinements)
global talkative

if talkative>1
    fprintf('\nStarting validation of simplices, time %s\n\n',datestr(now,13));
end

if nargin < 4 || isempty(bool_Hopf)
    bool_Hopf =0;
end
if nargin < 3 || isempty(index_validation_simplices)
    index_validation_simplices = length(list_of_simplices);
end
if nargin < 5 || isempty(max_refinements)
    max_refinements = 5;
end

ref = 0;
%while ref < max_refinements
    ref = ref + 1;
    [list_of_simplices, index_non_validated, Interval, Z0_iter, ...
        Z1_iter, Z2_iter, Y_iter] = iteration_validation(list_of_simplices,...
        list_of_nodes, index_validation_simplices, bool_Hopf);
%end

end

function [list_of_simplices, index_non_validated, Interval, Z0_iter, ...
    Z1_iter, Z2_iter, Y_iter] = iteration_validation(list_of_simplices,...
    list_of_nodes, index_validation_simplices, bool_Hopf)
global use_intlab

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
    use_intlab = 1;
    [flag,Imin,Imax,Yvector,Z0vector,Z1vector,Z2vector,...
        simplex_jj, list_of_nodes] = ...
        radii_polynomials_simplex(simplex_jj, list_of_nodes,bool_Hopf);
    use_intlab = 0;
    if flag < 1
        index_non_validated(end+1) = j;
    elseif talkative>0
        fprintf('The validation of the %i-th simplex succeeded\n', j)
    end
    if flag
        simplex_jj.verified = 1;
        list_of_simplices.simplex{j} = simplex_jj;
        %storage
        Interval(:,j)  = [Imin,Imax]';
        Z0_iter(:,j)   = vert(Z0vector);
        Z1_iter(:,j)   = vert(Z1vector);
        Z2_iter(:,j)   = vert(Z2vector);
        Y_iter(:,j)    = vert(Yvector);
    end
end
end