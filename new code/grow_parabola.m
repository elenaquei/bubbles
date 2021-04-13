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

% the connection matrix keep track of the edges
Connection_matrix = zeros(7,7);
Connection_matrix(1,:) = 1;
Connection_matrix(:,1) = 1;
Connection_matrix = Connection_matrix + diag(ones(6,1), -1) + diag(ones(6,1), 1);
Connection_matrix(2,end) = 1;
Connection_matrix(end,2) = 1;

% the simplex matrix keeps track of the filled simplices
Diag_and_over = diag(ones(6,1))+ diag(ones(5,1), 1);
Simplex_matrix = [ ones(6,1), Diag_and_over];
Simplex_matrix(end,2) = 1;

% second circle of points
for theta = pi/3 *([0:5]+0.5)
    y(end+1,:) = [sqrt(2)*sin(theta), sqrt(2)*cos(theta), 2];
    outer_points(end+1,:) = y(end,:);
end

plot_simplex(y, Connection_matrix, Simplex_matrix)

% growing 6 more 
for i = 1:0
    for k = 1:length(outer_points(:,1))
        arr = [outer_points(k,:);center];
        distArr(k) = pdist(arr,'euclidean');
    end
    [~,idx] = min(distArr);
    idx = idx(1); % if multiple, pick the first one
    % display the selected point on plot
end

% point = interractive_test(y);

function [y4] = predictor_corrector_alg(f, y_hat, y1, y2, y3)
F = @(y) [ f(y); dot(y- y_hat, y2-y1); dot(y-y_hat, y3-y1)]; 
% F is square - Newton works

DF =@(y) fin_diff(F, y, length(y1));

y4 = Newton(F, DF, y_hat);

end

function Df = fin_diff(f, x, length_x)

Df = zeros(length_x, length_x);
h_scal = 10^-5;

Id = eye(length_x);
f_x = f(x);

for i = 1: length_x
    h = h_scal * Id(:,i);
    f_x_plus_h = f(x + h);
    
    for j = 1: length_x
        Df(j,i) = ( f_x_plus_h(j) - f_x(j)) / h;
    end
    
end

end

function plot_simplex(point, connection_mat, simplex_mat)
figure()
plot(point(:,1), point(:,2),'o');
hold on
% plot the connections
for i=1:size(connection_mat,1)
    pointi = point(i,:);
    other_points = find(connection_mat(i,:));
    pointj = point(other_points,:);
    plot_lines(pointi, pointj)
end

% plot the simplices
for i = 1:length(simplex_mat(:,1))
    points = point(find(simplex_mat(i,:)),:);
    fill(points(:,1),points(:,2),'g')
end
end

function plot_lines(x,points)
for j = 1:size(points,1)
    y = points(j,:);
    plot([x(1), y(1)], [x(2), y(2)]);
end
end