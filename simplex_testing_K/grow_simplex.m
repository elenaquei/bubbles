function [x_new,R2frame,yframe] = grow_simplex(xc,xc_front,xc_int,xc_h)
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
gap_min = pi/6; % This is the minmum in the paper of G/L/P.
tau = 1.2;      % This is the tau used in the paper of G/L/P.
div_tol = 1E-2; % Expansion/divergence criteria for step size reduction.
tol = 1E-12;    % Tolerance for convergence of Gauss-Newton.
F = @(X) X(1) - X(2)^2 - X(3)^2;
DF = @(X) [1, -2*X(2), -2*X(3)];
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
% Get gap angle in local coordinate system and infer orientation. POSSIBLE
% BUG HERE.
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
    return
elseif n_new_simplex == 1   % New simplex is formed by extant nodes.
    x_new = [xc,xc_front];
    R2frame = [y1_R2,y0_R2,y2_R2];
    yframe = [yf(:,1),y0,yf(:,2)];
else
    % Generate predictor "fan" in R2 coordinate system. %%POSSIBLE BUG
    % HERE.
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