function [DxFn,Dalphafn,DxxFnV,DalphaxFn,DalphaxxFn,DalphaalphaFn,...
    DalphaalphaxFn,DxxxFnV] = derivatives_Hamiltonian(x,alpha,v)
%function [DxFn,Dalphafn,DxxFnV,DalphaxFn,DalphaxxFn,DalphaalphaFn,...
%    DalphaalphaxFn,DxxxFnV] = derivatives_Hamiltonian(x,alpha,v)


mu = alpha;
a = 2;
b = 3;
c = 1;
d = -1;

X = x(1); Y = x(2); Z = x(3); U = x(4);

DxFn =[ 0               1       0        0
        0               0       1        0
        0               0       0        1
  b+2*c*X+3*d*X^2       mu      a        0];

Dalphafn = [0
    0
    0
    Y];

if nargout<3
    return
end

v1=v(1); v2=v(2); v3=v(3); v4 = v(4);

DxxFnV = [0 0 0 0 
    0 0 0 0 
    0 0 0 0 
    4*c*v1+6*d*X*v1   0  0  0];

DalphaxFn = [0 0 0 0 
    0 0 0 0 
    0 0 0 0 
    0   1  0  0];

DalphaxxFn=zeros(4);
DalphaalphaFn=zeros(4,1);
DalphaalphaxFn=zeros(4);

DxxxFnV=[0 0 0 0 
    0 0 0 0 
    0 0 0 0 
    6*d*v1^2   0  0  0];