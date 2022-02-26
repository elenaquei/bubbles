function [big_X1,G_big, const_G_big,Es1,const_Es1_out] = archlength_saddle(X, G, const_G, Es0, const_Es0,Es1,const_Es1, step_size)
% function [X1,Es1] = archlength_saddle(X, G, const_G, Es, const_Es, h)
%
% predicts new point based on previous one
% INPUT
% X     previous point Hopf
% G     scalar equations (2)
% const_G   constants of the scalar equations
% Es    continuation equation
% const_Es  constant of the continuation equation
% h     step size
% OUTPUT
% X1    new point (after Newton)
% Es1   new coninuation equation
% const_Es1     constant of the new continuation equation


% Es = [Es; 0*Es; 0*Es];
% EDelta = [EDelta;0*EDelta;0*EDelta];
% x_tilde = PDE_vec2vec(X) + Es.'*h;
% 
% const_Es = x_tilde*Es;
% 
% x_tilde = PDE_vector(x_tilde, X.size_real, X.size_vector, X.node_space);
% 
% X1 = Newton_saddle(x_tilde,G,const_G,Es,const_Es, EDelta, const_EsDelta);
% %X1 = x_tilde;
% 
% [x1, y1, z1 ] = split(X1);
% 
% [Es1,const_Es1] = cont_equation(x1, G);
% 
% X1.parameters = X.parameters + rand(size(X.parameters)); % REMOVE
% 
% %[Es1,const_Es1] = cont_equation_saddle(X1,G,Es,EDelta);

[x1, ~, ~] = split(X);

[X1,Es1,const_Es1] = archlength(x1, G, const_G, Es1, const_Es1, step_size);
big_X1 = compute_saddle_from_Hopf(X1,G,Es0,Es1,const_Es0, const_Es1);
[G_big, const_G_big] = saddle_equations(G,const_G,Es0,const_Es0,Es1,const_Es1);
if size(Es1,1)>1
    Es1 = [Es1; 0*Es1; 0*Es1]; % continuation equation
else
    Es1 = [Es1, 0*Es1, 0*Es1];
end
const_Es1_out = const_Es1;
