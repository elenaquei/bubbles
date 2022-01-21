function merge_validated_results(file1, file2)
% function merge_validated_results(file1, file2)
%
% we assume the problems are the same, just two validation runs out of the
% same non-validte data
% save intersting info in file1

load(file1)

list_of_simplices1 = list_of_simplices;
step_size1 = step_size;
bool_Hopf1 = bool_Hopf;
list_of_nodes1 = list_of_nodes;
list_of_frontal_nodes1 = list_of_frontal_nodes;
bool_validated1 = bool_validated;
if bool_validated1
    Interval1 = Interval;
    Z0_iter1 = Z0_iter;
    Z1_iter1 = Z1_iter;
    Z2_iter1 = Z2_iter;
    Y_iter1 = Y_iter;
end

load(file2);
list_of_simplices2 = list_of_simplices;
step_size2 = step_size;
bool_Hopf2 = bool_Hopf;
list_of_nodes2 = list_of_nodes;
list_of_frontal_nodes12 = list_of_frontal_nodes;
bool_validated2 = bool_validated;
if bool_validated2
    Interval2 = Interval;
    Z0_iter2 = Z0_iter;
    Z1_iter2 = Z1_iter;
    Z2_iter2 = Z2_iter;
    Y_iter2 = Y_iter;
end


% first basic check
if step_size1~= step_size2 || bool_Hopf1~=bool_Hopf2 || any(list_of_frontal_nodes1~=list_of_frontal_nodes2)
    error('These files refer to different manifolds')
end
if length(list_of_simplices1)~=length(list_of_simplices2)
    error('Not the same number of simplices')
end

% if trivial cases, simply save and return
if (bool_validated2 && ~bool_validated1) ||(~bool_validated2 && ~bool_validated1)
    save(save_file,'list_of_simplices','list_of_nodes','Interval','Z0_iter',...
            'Z1_iter','Z2_iter','Y_iter','step_size','bool_Hopf', 'bool_validated',...
            'list_of_frontal_nodes');
    return
elseif ~bool_validated2 && bool_validated1
    return
end

% store all interesting data on file1
for i = 1:length(list_of_simplices1)
    simplex1 = list_of_simplices1.simplex{i};
    simplex2 = list_of_simplices2.simplex{i};
    if simplex1.verified && simplex2.verified
        Interval1(i,1) = min(Interval1(i,1), Interval2(i,1)); 
        Interval1(i,2) = max(Interval1(i,2), Interval2(i,2));
        simplex1.validation_interval = Interval1(i,:);
        Z0_iter1(i,:) = min(Z0_iter1(i,:),Z0_iter2(i,:));
        Z1_iter1(i,:) = min(Z1_iter1(i,:),Z1_iter2(i,:));
        Z2_iter1(i,:) = min(Z2_iter1(i,:),Z2_iter2(i,:));
        Y_iter1(i,:) = min(Y_iter1(i,:),Y_iter2(i,:));
    elseif ~simplex1.verified && simplex2.verified
        temp_simplex = simplex1;
        simplex1 = simplex2;
        simplex2 = temp_simplex;
        Interval1(i,:) = Interval2(i,:);
        Z0_iter1(i,:) = Z0_iter2(i,:);
        Z1_iter1(i,:) = Z1_iter2(i,:);
        Z2_iter1(i,:) = Z2_iter2(i,:);
        Y_iter1(i,:) = Y_iter2(i,:);
    end
    if simplex1.verified && ~simplex2.verified
        % do nothing
    elseif ~simplex1.verified && ~simplex2.verified
        % do nothing
    end
    list_of_simplices1.simplex{i} = simplex1;
end

list_of_simplices = list_of_simplices1;
bool_validated = bool_validated1;
if bool_validated1
    Interval = Interval1;
    Z0_iter = Z0_iter1;
    Z1_iter = Z1_iter1;
    Z2_iter = Z2_iter1;
    Y_iter = Y_iter1;
end

save(file1,'list_of_simplices','list_of_nodes','Interval','Z0_iter',...
            'Z1_iter','Z2_iter','Y_iter','step_size','bool_Hopf', 'bool_validated',...
            'list_of_frontal_nodes');