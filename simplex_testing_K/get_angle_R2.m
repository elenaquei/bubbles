function theta = get_angle_R2(X)
% INPUT X, vector in R2
% OUTPUT theta, positve scalar clockwise angle from axis e1.
x = X(1);   y = X(2);
theta = -atan2(y,x);
if x>=0
    if y<=0
        %theta = theta; % do nothing
    else
        theta = theta + 2*pi;
    end
elseif x<0
    if y<0
        %theta = theta; % do nothing
    else
        theta = theta + 2*pi;
    end 
end

