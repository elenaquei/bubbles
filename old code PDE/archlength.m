function [X1,Es,const_Es1, G1] = archlength(X, G, const_G, Es, const_Es0, h)
% function [X1,Es1] = archlength(X, G, const_G, Es, const_Es, h)
%
% predicts new point based on previous one
% INPUT
% X     previous point
% G     scalar equations (2)
% const_G   constants of the scalar equations
% Es    continuation equation
% const_Es  constant of the continuation equation
% h     step size
% OUTPUT
% X1    new point (after Newton)
% Es1   new continuation equation
% const_Es1     constant of the new continuation equation

x_tilde = PDE_vec2vec(X) + Es.'*h;

const_Es1 = const_Es0 - h;%- x_tilde*Es;

x_tilde = PDE_vector(x_tilde, X.size_real, X.size_vector, X.node_space);

%[Es1,const_Es1] = cont_equation(x_tilde,G);

X1 = Newton_HopfKS(x_tilde,G,const_G,Es,const_Es1);

%const_Es1 = X1.parameters(3);

% with the new continuation equation, only the constant is changed
G1 = new_amplitude(X1);

X1 = Newton_HopfKS(x_tilde,G1,const_G,Es,const_Es1);

end