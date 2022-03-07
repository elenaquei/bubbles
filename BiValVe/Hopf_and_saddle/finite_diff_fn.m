function [Dfn]=finite_diff_fn(fn,x,alpha)

h=1e-5;
M=length(x);
E=eye(M);
Dfn=zeros(M);
for j=1:M
    xh=x+h*E(:,j);
    fnxh=feval(fn,xh,alpha); fnx=feval(fn,x,alpha);
    Dfn(:,j)=(fnxh-fnx)/h;
end
end