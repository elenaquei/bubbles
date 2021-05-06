%function work_with_simplices()
%
% goal: grow a manifold of solutions of f=0 in 3D
%
% constraints: start from 7 points chosen by hand, with a "center", grow
% according to "small gap angle first".

% the zero finding problem
f=@(a,b,x) x - a^2 - b^2;

% or
fy = @(y) y(3) - y(1)^2 - y(2)^2;
DF =@(y) finite_diff(fy, y);

project([1;2;3],fy,DF);

center = [0,0,0];

y = center;
outer_points = [];

% first circle of points
for theta = pi/3 *[0:5]
    y(end+1,:) = [sin(theta), cos(theta), 1];
    outer_points(end+1,:) = y(end,:);
end

for i = 0:10
    theta = pi/11*i;
    angle_loc = angle([cos(theta),sin(theta)]);
    if abs(theta - angle_loc)>10^-6
        error('nope nope - angle wrong!')
    end
end

list_of_nodes = y;
list_of_simplices = cell(6,1);
simplex_1 = create_simplex(1,0,1, [1,2,3]);
simplex_2 = create_simplex(2,0,1, [1,3,4]);
simplex_3 = create_simplex(3,0,1, [1,5,4]);

for i = 1:5
    list_of_simplices{i} = create_simplex(i,0,1, [1,1+i,2+i]);
end
list_of_simplices{6} = create_simplex(6,0,1, [1,7,2]);


list_of_nodes = cell(size(y,1),1);
for i = 1:size(y,1)
    list_of_nodes{i} = create_node(i, y(i,:), i>=1);
end

list_of_nodes = update_patches(list_of_nodes,list_of_simplices);

% plot_simplex(simplex_1,list_of_nodes)
% hold on
% plot_simplex(simplex_2,list_of_nodes)
% plot_simplex(simplex_3,list_of_nodes)

figure;
plot_list_of_simplices(list_of_simplices,list_of_nodes)
hold on
plot_patch(list_of_nodes{3},list_of_simplices,list_of_nodes)

list_of_frontal_nodes = 2:7;
for i = 1:7
    if any(list_of_frontal_nodes==i) && list_of_nodes{i}.frontal~=1
        error('1');
    end
end
for i = 1:6
    if list_of_nodes{list_of_frontal_nodes(i)}.frontal~=1
        error('2')
    end
end

% figure;
% plot_list_of_simplices(list_of_simplices,list_of_nodes)
% hold on
% plot_patch(list_of_nodes{6},list_of_simplices,list_of_nodes)
% figure;
% plot_list_of_simplices(list_of_simplices,list_of_nodes)
% hold on
% plot_patch(list_of_nodes{7},list_of_simplices,list_of_nodes)

close all
[k,list_of_nodes,list_of_simplices,list_of_frontal_nodes] =...
    equate(6,7,list_of_nodes,list_of_simplices,list_of_frontal_nodes);
plot_list_of_simplices(list_of_simplices,list_of_nodes)

% angle( [9,0])
% angle ([1,1]) - pi/4
% angle( [0,-1]) - 3*pi/2 

coords_5 = project([0,-0.5,0.5]',fy,DF);
list_of_nodes{5}.coordinates = coords_5';

plot_list_of_simplices(list_of_simplices,list_of_nodes)

% growing the manifold
gap_angle_vec = 0 * list_of_frontal_nodes;
clockwise = gap_angle_vec;
for i=1:length(list_of_frontal_nodes)
    node_index = list_of_frontal_nodes(i);
    node = list_of_nodes{node_index};
    DF_x = DF(node.coordinates);
    [gap_angle_vec(i),clockwise(i)] = gap_angle(node, DF_x, list_of_simplices, list_of_nodes);
end
disp(list_of_frontal_nodes)
disp(gap_angle_vec)

[min_gap, priority_node_index] = min(gap_angle_vec);
priority_node = list_of_frontal_nodes(priority_node_index);
clockwise_ = clockwise(priority_node_index);

[list_of_nodes, list_of_simplices] = grow_patch(list_of_nodes{priority_node},...
    min_gap, clockwise_, list_of_nodes, list_of_simplices, fy, DF);



% %
% STRUCTURE
% %
% TODO list:
% - use gap angle to determine priority
% - use gap angle to propose new approximations
% - add the new nodes appropriately in new simplices

% works 
function x = project(x_approx, f, df)
% takes an approximation x_approx and projects it on the manifold f(x)=0

tol = 10^-6;
DF_x = df(x_approx);

[Q1Q2,R0] = qr(DF_x.');

Q1 = Q1Q2(:,1:end-2);
% Q2 = Q1Q2(:,end-2:end);
R = R0(1:end-2,:);

x = x_approx;

for iter = 1:100
    z = R\f(x);
    step = Q1 * z;
    x = x - step;
    if norm(step)<tol
        break
    end
end
if norm(step)>tol
    warning('Did not converge')
end

end


% works
function node = create_node(number, coords, frontal)
node.coordinates = coords;
node.number = number;
node.patch = [];
node.frontal = frontal;
end


% works - a bit more cumbersome than expected! 
function [i_loc,list_of_nodes, list_of_simplices,list_of_frontal_nodes] = ...
    equate(i_loc,j_loc,list_of_nodes, ...
    list_of_simplices, list_of_frontal_nodes)
% replace nodes i and j with their average, and return the number of the
% new node

if i_loc>j_loc
    temp = i_loc;
    i_loc = j_loc;
    j_loc = temp;
end

node_i = list_of_nodes{i_loc};
node_j = list_of_nodes{j_loc};

coords = 1/2 * (node_i.coordinates + node_j.coordinates);

if node_i.frontal ~= node_j.frontal
    warning('Merging of a frontal node with a non-frontal node')
end
new_node = node_i;
new_node.coordinates = coords; % this way, it maintains the patch too!
new_node.patch = [];
list_of_nodes{i_loc} = new_node;

for bigger_than_j = j_loc:length(list_of_nodes)-1
    list_of_nodes{bigger_than_j} = list_of_nodes{bigger_than_j+1};
    % list_of_nodes{bigger_than_j}.number = list_of_nodes{bigger_than_j}.number-1;
end
list_of_nodes(end) = [];

removal_list = [];
for index = 1:length(list_of_simplices)
    % replace the node j by the node i whenever necessary
    nodes_in_simplex = list_of_simplices{index}.nodes;
    index_node = nodes_in_simplex == j_loc;
    nodes_in_simplex(index_node) = i_loc;
    
    % target simplices with twice the same node
    if length(find(nodes_in_simplex==i_loc))>1
        removal_list(end+1) = index;
    end
    
    % % replace higher valued nodes by their new numbering
    % index_node = nodes_in_simplex >= j;
    % nodes_in_simplex(index_node) = nodes_in_simplex(index_node)-1;
    
    list_of_simplices{index}.nodes = nodes_in_simplex;
end

for index = removal_list
    for index_j = index:length(list_of_simplices)-1
        list_of_simplices{index_j} = list_of_simplices{index_j+1};
        list_of_simplices{index_j}.number = list_of_simplices{index_j}.number-1;
    end
    list_of_simplices(end) = [];
end

list_of_nodes = update_patches(list_of_nodes,list_of_simplices);

if nargin>4
    list_of_frontal_nodes = setdiff(list_of_frontal_nodes,j_loc);
end

end


function [list_of_nodes, list_of_simplices,list_of_frontal_nodes] = ...
    grow_patch(node, min_gap, clockwise, list_of_nodes, list_of_simplices, ...
    fy, DF,list_of_frontal_nodes)
% select open edges
out_nodes = open_edges(node, list_of_simplices);
if clockwise
    out_nodes = out_nodes(2:-1:1);% flip order
end


% decide what to do depending on min_gap
%       equate open edges
% or 
%       create new simplex
% or 
%       create new nodes
if min_gap<pi/6
    [~,list_of_nodes, list_of_simplices,list_of_frontal_nodes] = ...
    equate(out_nodes(1),out_nodes(2),list_of_nodes, ...
    list_of_simplices,list_of_frontal_nodes);
    return
   
elseif min_gap<pi/3
    % create new simplex
    number_simplex = length(list_of_simplices)+1;
    verified = (1==0);
    frontal = (1==1);
    list_of_simplices{end+1} = create_simplex(number_simplex,verified,...
        frontal, [node,out_nodes]); % DO we care about ordering?
else
    % create new points
    number_of_new_points = ceil(min_gap/pi*3);
    gap_angle = min_gap/number_of_new_points;
    for i = 1:number_of_new_points
        c
    end
end



end


% works
function list_of_nodes = update_patches(list_of_nodes,list_of_simplices)
if length(list_of_nodes)==1
    list_of_nodes_cell = cell(1,1);
    list_of_nodes_cell{1} = list_of_nodes;
    list_of_nodes = list_of_nodes_cell;
end
for j = 1: length(list_of_nodes)
    list_of_nodes{j}.patch = [];
    for k = 1: length(list_of_simplices)
        list_of_nodes{j} = add_simplex_to_patch(list_of_nodes{j},...
            list_of_simplices{k});
    end
end
if length(list_of_nodes)==1
    list_of_nodes = list_of_nodes{1};
end
end


% works
function node = add_simplex_to_patch(node, simplex)
if any(simplex.nodes == node.number)
    node.patch(end+1)= simplex.number;
end
end


% works
function simplex = create_simplex(number,verified,frontal, nodes)
simplex.number = number;
simplex.verified = verified;
simplex.fontal = frontal;
simplex.nodes = nodes;
simplex.verification_coeff = 0; % for later, "how easy it was to verify"
end


% works
function plot_simplex(simplex, list_of_nodes,color)
if nargin<3
    color = 'r';
end
x_coord = zeros(1,3);
y_coord = zeros(1,3);
z_coord = zeros(1,3);
label = cell(3,1);
for i = 1:3
    x_coord(i) = list_of_nodes{simplex.nodes(i)}.coordinates(1);
    y_coord(i) = list_of_nodes{simplex.nodes(i)}.coordinates(2);
    z_coord(i) = list_of_nodes{simplex.nodes(i)}.coordinates(3);
    label{i} = 'x'+string(list_of_nodes{simplex.nodes(i)}.number);
end
fill3(x_coord,y_coord,z_coord,color)
alpha 0.5
text(x_coord,y_coord,z_coord,label,'VerticalAlignment','bottom','HorizontalAlignment','right')
end

% works
function plot_list_of_simplices(list_simplex, list_of_nodes)
for i = 1: length(list_simplex)
    plot_simplex(list_simplex{i},list_of_nodes);
    hold on
end
hold off
end

% works
function plot_patch(node,list_simplex,list_of_nodes)
for i = 1: length(node.patch)
    plot_simplex(list_simplex{node.patch(i)},list_of_nodes,'b');
    hold on
end
hold off
end

% works = from simplex_1994
function alpha = angle_old(vector)
if length(vector)~=2
    error('Needs a projection on R2')
end
if norm(vector) ==0
    alpha = 0;
    return
end
a = vector(1); b = vector(2);
if b==0
    phi =0;
else
    phi = atan(b/a);
end
if abs(a) > abs(b)
    if a >=0
        alpha = phi;
        if b<0
            alpha = 2*pi-phi;
        end
    else
        alpha = pi - phi;
        if b<0
            alpha = pi + phi;
        end
    end
else
    if a==0
        phi =0;
    else
        phi = atan(a/b);
        if a<0 
            phi = -phi;
        end
    end
    if b>=0
        alpha = pi/2 - phi;
    else
        alpha = 3*pi/2 + phi;
    end
end
end

% works - tested more on R^2 
function alpha = angle(vector)

if norm(vector) ==0
    alpha = 0;
    return
else
    vector = vector/norm(vector);
end

if any(vector<0)
    if vector(1)<0 && vector(2)>0
        rotate_angle = pi/2;
    elseif vector(1)<0 && vector(2)<0
        rotate_angle = pi;
    else
        rotate_angle = 3*pi/2;
    end
else
    rotate_angle = 0;
end

rotation_matrix = [ cos(rotate_angle), sin(rotate_angle);
    -sin(rotate_angle), cos(rotate_angle)];

new_vector = rotation_matrix*[vector(1);vector(2)];

if new_vector(2) == 0
    alpha = 0+rotate_angle;
    return
end
alpha = abs(atan(new_vector(2)/new_vector(1)))+rotate_angle;

end



function [out_nodes,in_nodes] = open_edges(node, list_of_simplices)
out_nodes = [];
in_nodes = [];
for i = node.patch
    for j = list_of_simplices{i}.nodes
        if j == node.number
            continue
        end
        if any(out_nodes == j)
            out_nodes(out_nodes==j) = []; % if it's in, cancel it
        else
            out_nodes(end+1) = j;
        end
        in_nodes(end+1) = j;
    end
end
end

% from siplex 1994
function [gap,clockwise] = gap_angle(node, DF_x, list_of_simplices, list_of_nodes)

U = null(DF_x);

[out_nodes,in_nodes] = open_edges(node, list_of_simplices);

if length(out_nodes)~=2
    error('there are not two outnodes?!')
end
internal_nodes = setdiff(in_nodes,out_nodes);
if isempty(internal_nodes)
    internal_coord = 1/2 *(list_of_nodes{out_nodes(1)}.coordinates + ...
        list_of_nodes{out_nodes(2)}.coordinates); % average of two open edges
else
    internal_coord = list_of_nodes{internal_nodes(1)}.coordinates-node.coordinates;
end
internal_coord = (internal_coord) / norm(internal_coord); 
% open edges
edge1_in_Rn = list_of_nodes{out_nodes(1)}.coordinates;
edge1_in_Rn = (edge1_in_Rn - node.coordinates) / norm(edge1_in_Rn);
edge2_in_Rn = list_of_nodes{out_nodes(2)}.coordinates;
edge2_in_Rn = (edge2_in_Rn - node.coordinates) / norm(edge2_in_Rn);

if size(edge1_in_Rn,1) == 1
    edge1_in_Rn=edge1_in_Rn.';
end
if size(edge2_in_Rn,1) == 1
    edge2_in_Rn=edge2_in_Rn.';
end
if size(internal_coord,1) == 1
    internal_coord=internal_coord.';
end

edge1 = U.' * edge1_in_Rn;
edge2 = U.' * edge2_in_Rn;
internal_dir = U.' * internal_coord;

rotation_angle = angle(edge2);
rotation_matrix = [ cos(rotation_angle), sin(rotation_angle);
    -sin(rotation_angle), cos(rotation_angle)];

rotated_edge1 = rotation_matrix * edge1;
rotated_internal_dir = rotation_matrix * internal_dir;

angle_edge1 = angle(rotated_edge1);
angle_internal_dir = angle(rotated_internal_dir);

if angle_edge1 > angle_internal_dir
    gap = 2 * pi - angle_edge1;
    clockwise = (1==1);
else
    gap = angle_edge1;
    clockwise = (1==0);
end

% figure
% plot([0,1],[0,0])
% hold on
% plot([0,rotated_edge1(1)],[0,rotated_edge1(2)],'r')
% plot([0,rotated_internal_dir(1)],[0,rotated_internal_dir(2)],'g')
% title(string(node.number)+' angle:'+string(gap))

end