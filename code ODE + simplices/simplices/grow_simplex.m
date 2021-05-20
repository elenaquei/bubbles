function [list_of_nodes,list_of_simplices, list_of_new_frontal_nodes] = ...
    grow_simplex(node, step_size, list_of_nodes, list_of_simplices, problem)
% function [list_of_nodes,list_of_simplices] = grow_simplex(node, 
%           step_size, list_of_nodes, list_of_simplices, problem)
% 
% INPUT
% node          instance of the node class - the node that needs "growing"
% step_size     expected length of the new edges
% list_of_nodes cell - all old nodes
% list_of_simplices     instance of the all_simplices class - list of all
%               the simplices (at leat all the once in the patch of node)
% problem       instance of the problem class - the ODE we are solving
%
% OUTPUT
% list_of_nodes     cell - new list of nodes
% list_of_simplices    instance of the all_simplices class - list of all
%                   the simplices (old + new)



% old version: [x_new,R2frame,yframe] = grow_simplex(xc,xc_front,xc_int,xc_h)
%
% INPUT 
% xc : node to grow from.
% xc_front : frontal nodes incident to xc; there should be 2.
% xc_int : interior nodes incident to xc; should be at least 1.
% xc_h : length of shortest edge incident to xc in the database.
% OUTPUT
% x_new : matrix with columns xc followed by an ordered set of vertices such
% that x_new(:,2) and x_new(:,end) are from xc_front and all 'interior' 
% columns are computed by projecting these onto the tangent plane,
% computing a "fan" of new points between their gap angle, and projecting
% back onto the manifold.
% R2frame : The "fan", frontal nodes and gap complement direction after
% projection onto the tangent plane and subsequent (non-canonical)
% R2 (planar) isomorphism.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% wrapper
xc_h = step_size;
xc = Xi_vec2vec(node.solution);
[out_nodes, all_nodes] = open_edges(node, list_of_simplices);
in_nodes = setdiff(all_nodes, out_nodes);
if length(out_nodes) ~= 2
    error('There should be 2 and only 2 frontal edges')
end
xc_front = zeros(length(xc), length(out_nodes));
for i = 1:2
    xc_front(:,i) = Xi_vec2vec(list_of_nodes{out_nodes(i)}.solution);
end
xc_int = zeros(length(xc), length(in_nodes));
for i = 1:length(in_nodes)
    xc_int(:,i) = Xi_vec2vec(list_of_nodes{in_nodes(i)}.solution);
end

gap_min = pi/6; % This is the minmum in the paper of G/L/P.
tau = 1.2;      % This is the tau used in the paper of G/L/P.
div_tol = 1E-2; % Expansion/divergence criteria for step size reduction.
tol = 1E-12;    % Tolerance for convergence of Gauss-Newton.
F = @(X) Xi_vec2vec(apply(problem, vec2Xi_vec(X, node.solution))); 
DF = @(X) derivative_to_matrix(derivative(problem,vec2Xi_vec(X, node.solution),0)); 
[dim,n] = size(xc_int);

% Gap complement direction.
if n==1
    a = xc_int - xc;
else
    a = sum(xc_int,2)/n - xc;
    % a = xc_int(:,randi(n));   % Or pick a random interior one?
end
% Get projector
[Q,~] = qr(DF(xc).');
P = Q(:,end-1:end)*Q(:,end-1:end).';
% Project (xf-xc) for all fronal nodes xf, and complement avg., onto TM.
yf = P*(xc_front-xc);
yf(:,1) = yf(:,1)/norm(yf(:,1),2);  yf(:,2) = yf(:,2)/norm(yf(:,2),2);
y0 = P*a;   y0 = y0/norm(y0,2);
% Build R2 coordinate system around yf(:,1) and y0.
R2b1 = y0;  R2b2 = yf(:,1) - dot(R2b1,yf(:,1))*R2b1;    
R2b2 = R2b2/norm(R2b2,2);
TM_to_R2 = [R2b1.'; R2b2.'];    % Tangent space to R2 transformation
R2_to_TM = [R2b1, R2b2];        % Inverse transformation
y0_R2 = TM_to_R2*y0;            % Representatives of y0, yf in R2
y1_R2 = TM_to_R2*yf(:,1);
y2_R2 = TM_to_R2*yf(:,2);
% Compute angles relative to e1
angle_y0 = get_angle_R2(y0_R2);
angle_y1 = get_angle_R2(y1_R2);
angle_y2 = get_angle_R2(y2_R2);
% Compute angle beta1 clockwise from y1 to y2.
if angle_y1>angle_y2
    beta1 = 2*pi - (angle_y1 - angle_y2);
else 
    beta1 = angle_y2 - angle_y1;
end
% Compute angle beta2 clockwise from y0 to y2.
if angle_y0>angle_y2
    beta2 = 2*pi - (angle_y0 - angle_y2);
else
    beta2 = angle_y2 - angle_y0;
end
% Get gap angle in local coordinate system and infer orientation. 
if beta1<beta2
    gap = beta1;
    yor = 1;
else
    gap = 2*pi - beta1;
    yor = 2;
end
% Decide how many simplices are to be added to fill the gap.
n_new_simplex = max(1,round(gap/(pi/3)));
if gap < gap_min            % Gap is too small; merge simplices.
    x_new = [];
    R2frame = [];
    yframe = [yf(:,1),y0,yf(:,2)];
    % merge these two nodes
    [index_merged_node,list_of_nodes, list_of_simplices] = ...
    equate(list_of_nodes, out_nodes(1),out_nodes(2), ...
    list_of_simplices);
    list_of_new_frontal_nodes = [index_merged_node, -out_nodes(2)];
    for i = node.patch
        list_of_simplices{i}.frontal = 0;
    end
    return
elseif n_new_simplex == 1   % New simplex is formed by extant nodes.
    x_new = [xc,xc_front];
    R2frame = [y1_R2,y0_R2,y2_R2];
    yframe = [yf(:,1),y0,yf(:,2)];
    % create new simplex 
    nodes_number = [node.number, out_edges];
    simplex_number = length(list_of_simplices)+1;
    verified = 0;
    frontal = 1;
    
    simplex_x = simplex(nodes_number, simplex_number, verified, frontal,...
    	list_of_nodes);
    % add it to the list
    list_of_simplices = append(list_of_simplices, simplex_x);
    for i = node.patch
        list_of_simplices{i}.frontal = 0;
    end
    list_of_new_frontal_nodes = [];
else
    % Generate predictor "fan" in R2 coordinate system. 
    y_fan_R2 = zeros(2,n_new_simplex-1);
    for j=1:n_new_simplex-1
        theta = gap*j/n_new_simplex;
        Rtheta = [cos(theta),sin(theta); -sin(theta),cos(theta)];
        if yor == 1
            y_fan_R2(:,j) = Rtheta*y1_R2;
        elseif yor == 2
            y_fan_R2(:,j) = Rtheta*y2_R2;
        end
    end
    R2frame = [y1_R2,y0_R2,y2_R2,y_fan_R2];
    % Push the "fan" back into the tangent space T_{xc}M, and scale.
    y_fan = R2_to_TM*y_fan_R2;
    x_fan = zeros(dim,n_new_simplex-1);
    for j=1:n_new_simplex-1
        y_fan(:,j) = y_fan(:,j)/norm(y_fan(:,j),2);
        x_fan(:,j) = xc + y_fan(:,j)*xc_h*tau;
    end
    % Refine predictors with Gauss-Newton
    for j=1:n_new_simplex-1
       x_init = x_fan(:,j);
       h_local = xc_h;
       delta = inf;
       while delta>tol*(1+norm(x_init,2))
           if norm(x_init - x_fan(:,j))>div_tol   % Diverging.
               x_fan(:,j) = xc + y_fan(:,j)/norm(y_fan(:,j),2)*h_local;
               x_init = x_fan(:,j);
               h_local = h_local/tau;
           end
           x_fan(:,j) = GN(x_init,@(z)F(z),@(z)DF(z));
           delta = norm(x_fan(:,j)-x_init,2);
           x_init = x_fan(:,j);
       end
    end
    % Output nodes list
    x_new = [xc,xc_front(:,yor),x_fan,xc_front(:,mod(yor,2)+1)];
    yframe = [yf(:,1),y0,yf(:,2),y_fan];
    % TO DO : transform x_new into a list of new simplices
    list_new_nodes = zeros(1, n_new_simplex - 1);
    for i = 1:n_new_simplex - 1
        number = length(list_of_nodes)+1;
        solution = vec2Xi_vec(x_fan(:,i), node.solution);
        new_node = node(number, solution, problem);
        list_of_nodes{end+1} = new_node;
        list_new_nodes(i) = new_node.number;
    end
    frontal = 1;
    verified = 0;
    for i = 1:n_new_simplex
        number = length(list_of_simplices);
        if i ==1
            number1 = out_edges(yor);
        else
            number1 = list_new_nodes(i-1);
        end
        if i == n_new_simplex
            number2 = out_edges(mod(yor,2)+1);
        else
            number2 = list_new_nodes(i);
        end
        nodes_number = [node.number, number1, number2];
        simplex_new = simplex(nodes_number, number, verified, frontal,list_of_nodes);
        list_of_simplices = append(list_of_simplices, simplex_new);
    end
    list_of_new_frontal_nodes = list_new_nodes;
    list_new_nodes = union(list_new_nodes, out_nodes);
    list_of_nodes = update_patches(list_of_nodes,list_of_simplices, list_new_nodes);
    for i = node.patch
        list_of_simplices{i}.frontal = 0;
    end
end


end

function x_new = GN(x,F,DF)
[Q,R] = qr(DF(x).');
[dim,~] = size(x);
R = R(1:dim-2,1:dim-2);
Q1 = Q(:,1:end-2);
z = R\F(x);
x_new = x - Q1*z;
end