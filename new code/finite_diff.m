function [Dfn]=finite_diff(fn,x,varargin)

if size(x,1)==1
    x = x.';
end

h=1e-5;
M=length(x);
N=length(fn(x));
E=eye(M);
Dfn=zeros(N,M);
for j=1:M
    xh=x+h*E(:,j);
    if any(size(x)~=size(xh))
        error('incompatible sizes')
    end
    fnxh=feval(fn,xh,varargin{:}); 
    fnx=feval(fn,x,varargin{:});
    Dfn(:,j)=(fnxh-fnx)/h;
end
end