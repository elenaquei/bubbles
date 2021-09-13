function DF = fdiff_fix_p_R0(X,tau,L,K)
ep = 1E-2;
M = size(X)*[1;0]-2;
DF = zeros(M,M);
I = eye(M+2);
DF(:,1) = (function_F(X + ep*I(:,1),X,tau,L,K) - function_F(X,X,tau,L,K))/ep;
for n=2:M
    DF(:,n) = (function_F(X + ep*I(:,n+2),X,tau,L,K) - function_F(X,X,tau,L,K))/ep;
end

end