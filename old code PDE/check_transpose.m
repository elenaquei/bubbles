function x_trans=check_transpose(sol,input)
% function x_trans=check_transpose(sol,input)
% 
% this function checks that sol and input are in the same format, i.e. that
% the most important dimension is the same for both inputs and, in case it
% isn't, conforms sol to input. It returns sol in case they are or sol.' in
% case they aren't.
%
% Example    check_transpose(ones(4,3), rand(7,2)) return ones(4,3)
% Example2   check_transpose(rand(2,7), ones(4,3)) return rand(7,2)
 
global dim_bigger_than_nodes

if length(size(sol))~=2 || length(size(input))~=2 
    warning('this function just works with input of dimension 2')
    x_trans=sol;
    return
end

if all(size(sol)==size(input))
    x_trans=sol;
    return
end
if size(sol,1) == size(input,2) && size(input,1)== size(sol,2) 
    x_trans=sol.';
    return
end 

if ~dim_bigger_than_nodes
    n_nodes=(max(size(input))-1)/2;
%     if size(input,1)>size(input,2)
%         if size(sol,2)>1 && mod(size(sol,2)-1,n_nodes)==0
%             x_trans=sol;
%             return
%         else
%             x_trans=sol.';
%             return
%         end
%     else
%         if size(sol,2)>1 && mod(size(sol,2)-1,n_nodes)==0
%             x_trans=sol.';
%             return
%         else
%             x_trans=sol;
%             return
%         end
%     end
    if size(sol,1) > size(sol,2) && size(input,1)>size(input,2)
        x_trans=sol;
        return
    elseif size(sol,1) < size(sol,2) && size(input,1)<size(input,2)
        x_trans=sol;
        return
    else
        x_trans=sol.';
    end
else
    n_nodes=(min(size(input))-1)/2;
    
    if size(input,1)>size(input,2)
        if size(sol,2)>1 && mod(size(sol,2)-1,n_nodes)==0
            x_trans=sol;
            return
        else
            x_trans=sol.';
            return
        end
    else
        if size(sol,2)>1 && mod(size(sol,2)-1,n_nodes)==0
            x_trans=sol.';
            return
        else
            x_trans=sol;
            return
        end
    end
%     if size(sol,1) > size(sol,2) && size(input,1)<size(input,2)
%         x_trans=sol;
%         return
%     elseif size(sol,1) < size(sol,2) && size(input,1)<size(input,2)
%         x_trans=sol;
%         return
%     else
%         x_trans=sol.';
%     end
end