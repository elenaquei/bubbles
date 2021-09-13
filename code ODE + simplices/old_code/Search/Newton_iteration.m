function Xn = Newton_iteration(X,tau,L,K)
Xn = X;
DF = fdiff(Xn,tau,L,K);
F = function_F(Xn,Xn,tau,L,K);
while norm(F)>1E-15 & norm(F)<1E5
   Xn = Newton_step(Xn,F,DF);
   F = function_F(Xn,Xn,tau,L,K);
   norm(F)
end
Xn(1:7) = real(Xn(1:7));
end