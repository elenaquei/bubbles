function Norm=NuNorm(a,nu)
% function Norm=NuNorm(a,nu)
% a in l1nu
global use_intlab
local_intlab=0;
if ~isnumeric(a) && isintval(a)
    if use_intlab==1
        local_intlab=1;
    else
        a=sup(a);
    end
end
if ~isnumeric(nu) && isintval(nu)
    nu=inf(nu);
end

M=max(size(a));
if max(size(a))~=size(a,1)
    a=a.';
end
m=(M-1)/2;
K=[-m:m]';
if ~local_intlab
Norm=sum(abs(a)./(nu.^abs(K)));
else
    
Norm=sum(sup(abs(a)./(nu.^abs(K))));
end
return