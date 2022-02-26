function [G, const_G] = saddle_equations(G,const_G,Es0,const_Es0,Es1,const_Es1)
%
% INPUT
% Hopf continuation equation elements in s = 1
% Es0, const_Es0
% Es1, const_Es1
% OUTPUT
% elements for hte continuation equation in the saddle problem
% Es, const_Es
if size(Es1,1)>1
    Es1 = Es1.';
end
if size(Es0,1)>1
    Es0 = Es0.';
end
G = [G,0*G,0*G];
E1 = [Es1, 0*Es1, 0*Es1];
E2 = [0*Es1, Es1 - Es0, 0*Es1];
E3 = 2*[0*Es1, Es1 - Es0, 0*Es1];

const = [const_Es1, const_Es1-const_Es0, 0];

G = [G;E1;E2;E3];
const_G = [const_G,const];