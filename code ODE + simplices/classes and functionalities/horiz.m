function y = horiz(x)
% function y=horiz(x)
%
% returns a horizontal vector, either x or x.'

if min(size(x))~=1
    error('Works just with arrays')
end
if length(x) == size(x,1)
    y=x.';
else
    y=x;
end
return