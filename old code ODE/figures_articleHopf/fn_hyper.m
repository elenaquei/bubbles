function [fn,DxFn,Dalphafn] = fn_hyper(x,alpha)
% function [fn,DxFn,Dalphafn] = fn_hyper(x,alpha)
%
% according to the system in ??

a = alpha;
e = 2;
b=1;
c=b;
d=10;

x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4);

fn = [ a*x1 - a*x2 - x2*x3 + x4
    - b*x2 + x1*x3
    - c*x3 + d* x1 + x1*x2
    - e* x1 - e* x2];

if nargout<2 
    return
end

DxFn =[a    -a-x3    -x2   1
    x3      -b       x1    0 
    d+x2    x1       -c    0
    -e      -e       0     0];

Dalphafn = [x1-x2
    0
    0
    0];