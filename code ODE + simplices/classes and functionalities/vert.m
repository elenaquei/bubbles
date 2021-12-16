function y = vert(x)
% function y=vert(x)
%
% returns a vertical vector, either x or x.'

if min(size(x))~=1
    error('Works just with arrays')
end
if length(x) == size(x,2)
    y=x.';
else
    y=x;
end
return