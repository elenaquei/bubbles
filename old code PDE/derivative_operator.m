function K = derivative_operator(time_or_space, z)
% function K = derivative_operator(time_or_space, z)
%
% INPUT
% time_or_space      integer 1 or 2, indicating the dimension in which the
%                    derivative act
% z                  vector on which the output should act
% OUTPUT
% K                  matrix, such that Kz is the requested derivative
% REMARK: no imaginary part added

if time_or_space ~= 1 && time_or_space~=2
    error('Unknown derivative operator')
end

if isintval(z)
    is_intval_loc = 1;
else
    is_intval_loc = 0;
end

if min(size(z))==1 && length(size(z))==2
    if time_or_space~=2
        error('1D vector required to have a derivative in the second dimension')
    end
    nodes = (length(z) -1)/2;
    K = diag(-nodes:nodes);
    if is_intval_loc
        K = diag(intval(-nodes:nodes));
    end
    return
end

nodes = (size(z) -1)/2;

k = -nodes(time_or_space):nodes(time_or_space);
if is_intval_loc
    k = (intval(-nodes(time_or_space):nodes(time_or_space)));
end

if time_or_space == 1
    
    K = repmat(k,1,1+2*nodes(2));
    
else
    
    K = reshape(repmat(k,1+2*nodes(1),1),1,[]);
    
end

K = diag(K);