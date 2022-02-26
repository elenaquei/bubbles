function [DxFn,Dalphafn,DxxFnV,DalphaxFn,DalphaxxFn,DalphaalphaFn,...
    DalphaalphaxFn,DxxxFnV] = derivatives_hyper(x,alpha,v)
%function [DxFn,Dalphafn,DxxFnV,DalphaxFn,DalphaxxFn,DalphaalphaFn,...
%    DalphaalphaxFn,DxxxFnV] = derivatives_hyper(x,alpha,v)
%
% 
% according to the system in ??

a = alpha;
e = 2;
b=1;
c=b;
d=10;

x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4);

DxFn =[a    -a-x3    -x2   1
    x3      -b       x1    0 
    d+x2    x1       -c    0
    -e      -e       0     0];

Dalphafn = [x1-x2
    0
    0
    0];
if nargout<3
    return
end

v1=v(1); v2=v(2); v3=v(3); v4 = v(4);

DxxFnV = [0  -v3   -v2   0
    v3    0    v1   0 
    v2    v1   0    0
    0     0    0    0];

DalphaxFn = [1  -1   0    0
    0     0    0    0 
    0     0    0    0
    0     0    0    0];

DalphaxxFn=zeros(4);
DalphaalphaFn=zeros(4,1);
DalphaalphaxFn=zeros(4);

DxxxFnV=zeros(4);