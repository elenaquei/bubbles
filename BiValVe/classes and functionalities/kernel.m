function [x,iter] = kernel(A,x0)

% This subroutine computes a (unit) basis for the kernel of A. The 
% algorithm assumes that dim(ker(A))=1. The second input x0 is optional 
% and is assumed to be an initial guess (of length 1) for the kernel. 

if size(A,1)==size(A,2)
    error('Matrix needs to be rectangular, not square');
end

% Constants.                          
Nmax = 20;
tol = 1e-12;                        

% Compute initial guess if not supplied by the user. 
if(nargin==1)
    [Q,~] = qr(A'); % Needs some more comments...
    x0  = Q(:,end)/norm(Q(:,end));
end

% Define maps. 
F=@(x)([A*x; x0'*x-1]);
DF = [A;x0'];

% Initialization. 
x = x0; 
F0 = F(x); 
error_iter = norm(F0);
iter = 0; 

% Newton's method.
while (error_iter > tol && iter < Nmax)
    x = x - DF\F0; 
    F0 = F(x); 
    error_iter = norm(F0,inf);
    iter = iter + 1;
end



end

