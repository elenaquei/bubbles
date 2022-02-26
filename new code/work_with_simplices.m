%function work_with_simplices()
%
% goal: grow a manifold of solutions of f=0 in 3D
%
% constraints: start from 7 points chosen by hand, with a "center", grow
% according to "small gap angle first".
global debug_bool
debug_bool = 1;

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
n = 4;
for theta = pi/n *[0:n*2-1]
    theta_rand = theta*(0.9 + 0*rand);
    y_new =  10^-1 *(0.9 + 0*rand) * [sin(theta_rand), cos(theta_rand), 1];
    y(end+1,:) = project(y_new.', fy, DF);
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
list_of_simplices = cell(n*2,1);
simplex_1 = create_simplex(1,0,1, [1,2,3]);
simplex_2 = create_simplex(2,0,1, [1,3,4]);
simplex_3 = create_simplex(3,0,1, [1,5,4]);

for i = 1:n*2-1
    list_of_simplices{i} = create_simplex(i,0,1, [1,1+i,2+i]);
end
list_of_simplices{n*2} = create_simplex(n*2,0,1, [1,n*2+1,2]);


list_of_nodes = cell(size(y,1),1);
for i = 1:size(y,1)
    list_of_nodes{i} = create_node(i, y(i,:), i>=1);
end

list_of_nodes = update_patches(list_of_nodes,list_of_simplices);

% plot_simplex(simplex_1,list_of_nodes)
% hold on
% plot_simplex(simplex_2,list_of_nodes)
% plot_simplex(simplex_3,list_of_nodes)

% figure;
% plot_list_of_simplices(list_of_simplices,list_of_nodes)
% hold on
% plot_patch(list_of_nodes{3},list_of_simplices,list_of_nodes)

list_of_frontal_nodes = 2:n*2+1;
for i = 1:n*2+1
    if any(list_of_frontal_nodes==i) && list_of_nodes{i}.frontal~=1
        error('1');
    end
end
for i = 1:n*2-1
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

% close all
% [k,list_of_nodes,list_of_simplices,list_of_frontal_nodes] =...
%    equate(6,7,list_of_nodes,list_of_simplices,list_of_frontal_nodes);
% plot_list_of_simplices(list_of_simplices,list_of_nodes)

% angle( [9,0])
% angle ([1,1]) - pi/4
% angle( [0,-1]) - 3*pi/2

% coords_5 = project([0,-0.5,0.5]',fy,DF);
% list_of_nodes{5}.coordinates = coords_5';

% plot_list_of_simplices(list_of_simplices,list_of_nodes)

% growing the manifold
for k = 1:4
    gap_angle_vec = 0 * list_of_frontal_nodes;
    for i=1:length(list_of_frontal_nodes)
        node_index = list_of_frontal_nodes(i);
        node = list_of_nodes{node_index};
        DF_x = DF(node.coordinates);
        gap_angle_vec(i) = gap_angle(node, DF_x, list_of_simplices, list_of_nodes);
    end
    disp(list_of_frontal_nodes)
    disp(gap_angle_vec)
    
    [~, priority_node_index] = min(abs(gap_angle_vec));
    priority_node = list_of_frontal_nodes(priority_node_index);
    min_gap = gap_angle_vec(priority_node_index);
    
    %figure
    plot_list_of_simplices(list_of_simplices,list_of_nodes)
    pause(1)
    
    [list_of_nodes, list_of_simplices,list_of_frontal_nodes] =...
        grow_patch(list_of_nodes{priority_node},min_gap, ...
        list_of_nodes, list_of_simplices, fy, DF,list_of_frontal_nodes);
    
    plot_list_of_simplices(list_of_simplices,list_of_nodes)
    
end
% %
% STRUCTURE
% %
% TODO list:
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


function x = Newton(x_approx, f, df)
x = x_approx;
tol = 10^-6;
for iter = 1:100
    step = df(x) \ f(x);
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
    grow_patch(node, min_gap, list_of_nodes, list_of_simplices, ...
    f, df,list_of_frontal_nodes, step_size)
% select open edges
global debug_bool

[out_nodes, in_nodes] = open_edges(node, list_of_simplices);
if 1==0 && min_gap<0
    out_nodes = out_nodes(2:-1:1);% flip order
    min_gap = abs(min_gap);
end
if nargin > 7
    h = step_size;
else
    h = norm(node.coordinates - list_of_nodes{out_nodes(1)}.coordinates);
end

% decide what to do depending on min_gap
%       equate open edges
% or 
%       create new simplex
% or 
%       create new nodes

number_center_node = node.number;

if abs(min_gap)<pi/6
    [~,list_of_nodes, list_of_simplices,list_of_frontal_nodes] = ...
    equate(out_nodes(1),out_nodes(2),list_of_nodes, ...
    list_of_simplices,list_of_frontal_nodes);
   
elseif abs(min_gap)<pi/3
    % create new simplex
    number_simplex = length(list_of_simplices)+1;
    verified = (1==0);
    frontal = (1==1);
    list_of_simplices{end+1} = create_simplex(number_simplex,verified,...
        frontal, [node.number,out_nodes]); % DO we care about ordering?
else
    % create new points
    verified = (1==0);
    frontal = (1==1);
    number_of_new_points = ceil(abs(min_gap)/pi*3)-1;
    indices_new_points = zeros(number_of_new_points,1);
    
    number_old_node = out_nodes(1);
    y_old = list_of_nodes{number_old_node}.coordinates;
    y_center = node.coordinates;
    
    proj_on_2 = @(u,v) dot(u,v)/dot(v,v)*v;
    
    tang_plane = null(df(y_center));
    U = [list_of_nodes{out_nodes(1)}.coordinates-y_center;
        list_of_nodes{in_nodes(1)}.coordinates-y_center].';
    
    % out nodes projected on the tangent plane
    w1 = proj_on_2(U(:,1),tang_plane(:,1))+proj_on_2(U(:,1),tang_plane(:,2));
    w2 = proj_on_2(U(:,2),tang_plane(:,1))+proj_on_2(U(:,2),tang_plane(:,2));
    
    w1 = w1/norm(w1);
    w2 = w2 - dot(w1,w2)*w2;
    w2 = -w2 / norm(w2);
    
    gap_angle_loc = min_gap /(number_of_new_points+2);
    s = linspace(0,1,number_of_new_points+2);
    
    for i = 1:number_of_new_points
        x_new = cos((i+1)*gap_angle_loc)*w1 + sin((i+1)*gap_angle_loc)*w2;
        x_new = x_new/norm(x_new);
        y_new = y_center + h * x_new.';
        if debug_bool
            hold on
            plot3([y_center(1),y_new(1)],[y_center(2),y_new(2)],[y_center(3),y_new(3)], 'd')
        end
        coord = project(y_new.', f, df);%extend_project(y_new.', f, df,y_center+tang_plane(:,1).', y_center+tang_plane(:,2).');
        number_new_node = length(list_of_nodes)+1;
        list_of_nodes{end+1} = create_node(number_new_node, coord.', frontal);
        number_simplex = length(list_of_simplices)+1;
        nodes_in_new_simplex = [number_new_node, number_old_node, number_center_node];
        % DO WE CARE ABOUT ORDER????
        list_of_simplices{end+1} = create_simplex(number_simplex,...
            verified,frontal, nodes_in_new_simplex);
        
        list_of_frontal_nodes(end+1) = number_new_node;
        % update in for loop
        indices_new_points(i) = number_new_node;
        number_old_node = number_new_node;
        y_old = y_new;
    end
    % last new simplex - reconnect to out edges
    number_simplex = length(list_of_simplices)+1;
    nodes_in_new_simplex = [out_nodes(2), number_old_node, number_center_node];
    % DO WE CARE ABOUT ORDER????
    list_of_simplices{end+1} = create_simplex(number_simplex,...
        verified,frontal, nodes_in_new_simplex);
    for  i = 1:length(indices_new_points)
        index = indices_new_points(i);
        list_of_nodes{index} = update_patches(list_of_nodes{index},list_of_simplices);
    end
    list_of_nodes{out_nodes(1)} = update_patches(list_of_nodes{out_nodes(1)},list_of_simplices);
    list_of_nodes{out_nodes(2)} = update_patches(list_of_nodes{out_nodes(2)},list_of_simplices);
end
% if procedure successful, starting point is not frontal anymore
list_of_frontal_nodes = setdiff(list_of_frontal_nodes,node.number);

end


function coord = extend_project(x_new, f, df, x1, x2)
if size(x1,1)~=1
    x1=x1.';
end
if size(x2,1)~=1
    x2=x2.';
end
if size(x_new,1)~=1
    x_new=x_new.';
end
F = @(x) [ dot(x,x1-x_new) - dot(x_new, x1-x_new);
    dot(x,x2-x_new) - dot(x_new, x2-x_new);
    f(x)];
DF = @(x) [x1;
    x2;
    df(x)];
coord = Newton(x_new.', F, DF);
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

function alpha = angle_2(vector1,vector2)
rotation_angle = angle(vector2);
rotation_matrix = [ cos(rotation_angle), sin(rotation_angle);
    -sin(rotation_angle), cos(rotation_angle)];

rotated_edge1 = rotation_matrix * vector1;
alpha = angle(rotated_edge1);

end

function [out_nodes,all_nodes] = open_edges(node, list_of_simplices)
out_nodes = [];
all_nodes = [];
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
        all_nodes(end+1) = j;
    end
end
all_nodes = unique(all_nodes);
end

% from simplex 1994 - but negative gap angle when clockwise
function [gap] = gap_angle(node, DF_x, list_of_simplices, list_of_nodes)

U = null(DF_x);

[out_nodes,all_nodes] = open_edges(node, list_of_simplices);

if length(out_nodes)~=2
    error('there are not two outnodes?!')
end
internal_nodes = setdiff(all_nodes,out_nodes);
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
    gap = -(2 * pi - angle_edge1);
else
    gap = angle_edge1;
end

% figure
% plot([0,1],[0,0])
% hold on
% plot([0,rotated_edge1(1)],[0,rotated_edge1(2)],'r')
% plot([0,rotated_internal_dir(1)],[0,rotated_internal_dir(2)],'g')
% title(string(node.number)+' angle:'+string(gap))

end