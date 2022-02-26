function bool = isint(a)
% function bool=isint(a)
%
% returns 1 if input is an integer, 0 otherwise

bool =0;
if size(a,1)==1 && size(a,2) ==1 && mod(a,1) ==0
    bool =1;
end
return
end