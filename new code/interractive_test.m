function [coord] = interractive_test(points)
x = points(:,1);
y = points(:,2);
% z = points(:,3);

n_points = length(x);

if length(y)~=n_points %||length(z)~=n_points
    error('All coordinates must have the same length')
end

handle.x = x;
handle.y = y;
% handle.z = z;
% plot in 3D
handle.a = axes;
handle.p = plot(handle.x,handle.y,'ro');
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis');
% add callback when point on plot object 'handle.p' is selected
% 'click' is the callback function being called when user clicks a point on plot
handle.p.ButtonDownFcn = {@click,handle};
disp('Select a point - you only have 4 seconds')
pause(4)
% definition of click
function click(~,~,handle)
    % co-ordinates of the current selected point
    Pt = handle.a.CurrentPoint(2,1:2);
    distArr = zeros(1, n_points);
    % find point closest to selected point on the plot
    for k = 1:n_points
        arr = [handle.x(k) handle.y(k);Pt];
        distArr(k) = pdist(arr,'euclidean');
    end
    [~,idx] = min(distArr);
    % display the selected point on plot
    coord = [handle.x(idx) handle.y(idx)];
end
end