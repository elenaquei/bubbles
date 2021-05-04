function DHopf = derivative_Hopf(x_PDE, G, Es, flagFullSize, flagOp)
% function DHopf = derivative_Hopf(x_PDE, G, Es, flagFullSize, flagOp)
%
% INPUT 
% x_PDE     PDE_vector with 3 unknown parameters, one 1D vector and one
%           2D vector
% Es        PDE_vector or vector compatible with the previous PDE_vector,
%           continuation equation
% G         matrix of 2 lines, each line compatible with the
%           PDE_vector
% Es        (DEFAULT empty) elements for the continuation equation
% flagFullSize (DEFAULT 0) returns as big as possible
% flagOp    (DEFAULT 1) returns an overoperator
% OUTPUT
% DHopf     matrix, derivative of Kuramoto-shivashinki Hopf in x_PDE
% 

% WORKS - tested

if x_PDE.size_real ~= 3
    error('Input does not respect dimensional contraints - parameters')
end
if x_PDE.size_vector ~= 1
    error('Input does not respect dimensional contraints - vector')
end
if nargin<3 || isempty(Es)
    Es_in = 0;
    if isintval(G(1))
        is_intval_loc = 1;
    else
        is_intval_loc = 0;
    end
else
    Es_in = 1;
    if isintval(Es(1)) || isintval(x_PDE)
        is_intval_loc = 1;
    else
        is_intval_loc = 0;
    end
end
if nargin<4 || isempty(flagFullSize)
    flagFullSize=0;
end
if nargin<5 || isempty(flagOp)
    flagOp = 0;
end

if flagFullSize
    G_new = zeros(size(G,1), 3 + 2*2*x_PDE.node_space+1+(2*2*x_PDE.node_time+1)*(2*2*x_PDE.node_space+1));
    if is_intval_loc
        G_new = intval(G_new);
    end
    for i =1:size(G,1)
        G_new(i,:) = PDE_vec2vec(reshape(PDE_vector(G(i,:),3,1,...
            x_PDE.node_space),2*x_PDE.node_space, 2*x_PDE.node_time));
    end
    G=G_new;
    Es = PDE_vec2vec(reshape(PDE_vector(Es.',3,1,...
            x_PDE.node_space),2*x_PDE.node_space, 2*x_PDE.node_time));
    x_PDE = reshape(x_PDE,x_PDE.node_space*2,x_PDE.node_time*2);
end

y = x_PDE.vector1D;
z2D = squeeze(x_PDE.vector2D);
z1D = z2D(:);
y2D = embedding(y,z2D);
lambda = x_PDE.parameters(2);
nu = x_PDE.parameters(1);
a = x_PDE.parameters(3);

if ~Es_in
    
    DHopf = zeros(size(G,2)-1,size(G,2));
    if is_intval_loc
        DHopf = intval(DHopf);
    end
    DHopf(1:x_PDE.size_real-1,:) = G(1:2,:);
    
else
    
    DHopf = zeros(size(G,2),size(G,2));
    if is_intval_loc
        DHopf = intval(DHopf);
    end
    DHopf(1,:) = Es(:).';
    DHopf(2:x_PDE.size_real,:) = G(1:2,:);
end

K_1D = derivative_operator(2,y);
K_2D = derivative_operator(2,z2D);
J_2D = derivative_operator(1,z2D);

F1_y = - nu * K_1D.^4 + K_1D.^2 + 2 * 1i * K_1D *  conv_operator( y ,1 );

%F1_y(x_PDE.node_space+1,:) =1; 

F1_nu = - K_1D.^4 * y.';
F1_lambda = 0 * F1_nu;
F1_a = 0 * F1_lambda;

non_lin = ( 2 * Convo(z2D,y) + a * Convo(z2D,z2D));
F2_lambda = (- nu * K_2D.^4 + K_2D.^2) * z1D + 1i * K_2D * non_lin(:);

F2_nu = - lambda * K_2D.^4 * z1D; 

z_conv_z = Convo(z2D,z2D);
F2_a = lambda * 1i * K_2D * z_conv_z(:);

F2_y = 2 * lambda * 1i * K_2D * conv_operator(z2D,1); 

column_selection = x_PDE.node_time+1+(0:2*x_PDE.node_space)*(2*x_PDE.node_time+1);
F2_y = F2_y(:,column_selection);

F2_z = - 1i * J_2D - lambda * nu * K_2D.^4 + lambda * K_2D.^2 ...
+ lambda * 1i * K_2D * ( 2 * conv_operator(y2D,1) + 2 * a * conv_operator(z2D,1));

F1_z = zeros(size(F1_y,1),size(F2_z,2));

row_G4_base =(2* x_PDE.node_space + 1) +  ((2* x_PDE.node_space + 1) *(1+2* x_PDE.node_time)+1)/2;


if ~Es_in
    DHopf(x_PDE.size_real:end,:)= [ F1_nu, F1_lambda, F1_a, F1_y, F1_z
        F2_nu, F2_lambda, F2_a, F2_y, F2_z];
    row_G3 = 2 +  x_PDE.node_space + 1;
    row_G4 = 2 + row_G4_base;
    
else
    DHopf(1+x_PDE.size_real:end,:)= [ F1_nu, F1_lambda, F1_a, F1_y, F1_z
        F2_nu, F2_lambda, F2_a, F2_y, F2_z];
    row_G4 = 3 + row_G4_base;
    row_G3 = 3 +  x_PDE.node_space + 1;
end

DHopf(row_G3,:) = G(3,:);
DHopf(row_G4,:) = G(4,:);

if size(DHopf,1)>10^3
    DHopf=sparse(DHopf);
end

if size(DHopf,1)==size(DHopf,2) && flagOp
    DHopf = compose_overoperator(DHopf, x_PDE);
end