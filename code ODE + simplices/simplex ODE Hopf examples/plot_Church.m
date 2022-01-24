function plot_Church(list_simplex, list_of_nodes, index_simplices, varargin)
if nargin <3 || isempty(index_simplices)
    index_simplices = 1:length(list_simplex);
end

for index = 1: length(index_simplices)
    i = index_simplices(index);
    if length(index_simplices)>10 && ~list_simplex.simplex{i}.frontal
        if length(varargin)>1
            plot_simplex(list_simplex.simplex{i}, list_of_nodes, varargin{1},'LineStyle','none',varargin{2:end});
        elseif length(varargin)==1
            plot_simplex(list_simplex.simplex{i}, list_of_nodes, varargin{1},'LineStyle','none');
        else
            plot_simplex(list_simplex.simplex{i}, list_of_nodes, [], 'LineStyle','none');
        end
    else
        plot_simplex(list_simplex.simplex{i}, list_of_nodes, varargin{:});
    end
    hold on
end
alpha 0.5
set(gca,'FontSize',18)
hold off
xlabel('$\alpha$','Interpreter','Latex', 'FontSize', 20);
ylabel('$\beta$','Interpreter','Latex', 'FontSize', 20);
zlabel('$z_0$','Interpreter','Latex', 'FontSize', 20);
end

function plot_simplex(simplex, list_of_nodes, color, varargin)
% function plot_simplex(simplex, color)
if nargin<3 || isempty(color)
    if simplex.verified
        color = 'b';
    elseif simplex.frontal
        color = 'y';
    else
        color = 'r';
    end
end
alpha_coord = zeros(1,3);
beta_coord = zeros(1,3);
z0_coord = zeros(1,3);
label = cell(3,1);
modes = list_of_nodes{simplex.nodes_number(1)}.solution.nodes;
for i = 1:3
    % alpha vs beta vs z(0)
    alpha_coord(i) = list_of_nodes{simplex.nodes_number(i)}.solution.scalar(2);
    beta_coord(i) = list_of_nodes{simplex.nodes_number(i)}.solution.scalar(3);
    z0_coord(i) = list_of_nodes{simplex.nodes_number(i)}.solution.vector(3,modes+1);
    label{i} = 'x'+string(simplex.nodes_number(i));
end
% label_simplex = 'S' +string(simplex.number);
% center_x = mean(p_coord);
% center_y = mean(R0_coord);
% center_z = mean(a_coord);
fill3(alpha_coord,beta_coord,z0_coord,color, varargin{:})

%text(x_coord,y_coord,z_coord,label,'VerticalAlignment','bottom','HorizontalAlignment','right')
%text(center_x,center_y,center_z,label_simplex,'VerticalAlignment','bottom','HorizontalAlignment','right')
end
