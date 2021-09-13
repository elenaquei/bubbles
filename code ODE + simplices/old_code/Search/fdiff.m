function DF = fdiff(X,tau,L,K)
ep = 1E-5;
M = size(X)*[1;0]-2;
DF = zeros(M,M+2);
I = eye(M+2);
for n=1:M+2
    DF(:,n) = (function_F(X + ep*I(:,n),X,tau,L,K) - function_F(X,X,tau,L,K))/ep;
end

end