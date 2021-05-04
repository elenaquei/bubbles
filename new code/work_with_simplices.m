%function work_with_simplices()
%
% goal: grow a manifold of solutions of f=0 in 3D
%
% constraints: start from 7 points chosen by hand, with a "center", grow
% accoring to "close to center algorithm". Every 3 steps, ask for new
% center interactively

% the zero finding problem
f=@(a,b,x) x - a^2 - b^2;

% or
fy = @(y) y(3) - y(1)^2 - y(2)^2;
DF =@(y) fin_diff(F, y, 3);

center = [0,0,0];

y = center;
outer_points = [];

% first circle of points
for theta = pi/3 *[0:5]
    y(end+1,:) = [sin(theta), cos(theta), 1];
    outer_points(end+1,:) = y(end,:);
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
% we good!


% angle( [9,0])
% angle ([1,1]) - pi/4
% angle( [0,-1]) - 3*pi/2 

% growing the manifold
for i=1:length(list_of_frontal_nodes)
    node = list_of_frontal_nodes{i};
    gap_angle_vec(i) = gap_angle(node, list_of_simplices, list_of_nodes);
end




% %
% STRUCTURE
% %
% TODO list:
% - compute gap angle
% - use gap angle to determine priority
% - use gap angle to propose new approximations
% - use Newton as corrector
% - add the new nodes appropriately in new simplices



% works
function node = create_node(number, coords, frontal)
node.coordinates = coords;
node.number = number;
node.patch = [];
node.frontal = frontal;
end

% works
function list_of_nodes = update_patches(list_of_nodes,list_of_simplices)
for j = 1: length(list_of_nodes)
    for k = 1: length(list_of_simplices)
        list_of_nodes{j} = add_simplex_to_patch(list_of_nodes{j},...
            list_of_simplices{k});
    end
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
for i = 1:3
    x_coord(i) = list_of_nodes{simplex.nodes(i)}.coordinates(1);
    y_coord(i) = list_of_nodes{simplex.nodes(i)}.coordinates(2);
    z_coord(i) = list_of_nodes{simplex.nodes(i)}.coordinates(3);
end
fill3(x_coord,y_coord,z_coord,color)
alpha 0.5
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
function alpha = angle(vector)
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


% from siplex 1994
function gap = gap_angle(node, list_of_simplices, list_of_nodes)
gap = 1;
end