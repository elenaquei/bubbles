function Z0_vector=Z0_bound(DHm,A,xBar)
% function Z0_vector=Z0_bound(DHm,A,xBar)
% 
% INPUT:
% DHm     complex matrix, can also be called Adagger, is the approximate 
%         derivative in the numerical solution xBar (given);
% A       complex matrix, numerical inverse of DHm;
% xBar    Xi_vector, numerical solution. Theoretically, there is no need
%         of having the numerical solution in this function, the goal is to have
%         the data if the problem inherent the solution (the size, the number
%         of nodes, ecc.)
%
% OUTPUT:
% Z0_vector   a standard vector of length = size_vector+size_scalar, the
%             values of this vector come from the Z0 bound in the associated pdf (we
%             refer to rigorous numerics for analytic solution od D.E..... for the
%             notation, mainly)

global talkative

B=eye(size(A))-A*DHm;

Z0_vector=norm_Ximat(B,xBar);


% error handling
if any(Z0_vector>1) %&& talkative>0
%    warning('Z1>1');
    fprintf('Z0 computed, %d\n',max(Z0_vector));
    %error('Z0 is bigger than 1, no interval found')
elseif talkative>2
   fprintf('Z0 computed, %d\n',max(Z0_vector));
   fprintf('\n');
elseif talkative>1
   fprintf('Z0 computed, %d\n',max(Z0_vector));
   fprintf('\n');
end
