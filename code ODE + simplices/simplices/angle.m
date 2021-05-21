function alpha = angle(vector)

if norm(vector) ==0
    alpha = 0;
    return
% else
%     vector = vector/norm(vector);
end

angle = atan2(vector(1),vector(2));

alpha = pi - angle;

% if any(vector<0)
%     if vector(1)<0 && vector(2)>0
%         rotate_angle = pi/2;
%     elseif vector(1)<0 && vector(2)<0
%         rotate_angle = pi;
%     else
%         rotate_angle = 3*pi/2;
%     end
% else
%     rotate_angle = 0;
% end
% 
% rotation_matrix = [ cos(rotate_angle), sin(rotate_angle);
%     -sin(rotate_angle), cos(rotate_angle)];
% 
% new_vector = rotation_matrix*[vector(1);vector(2)];
% 
% if new_vector(2) == 0
%     alpha = 0+rotate_angle;
%     return
% end
% alpha = abs(atan(new_vector(2)/new_vector(1)))+rotate_angle;

end