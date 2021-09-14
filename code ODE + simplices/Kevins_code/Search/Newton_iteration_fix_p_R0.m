function Xn = Newton_iteration_fix_p_R0(X,tau,L,K,N)
Xn = X;
DF = function_DF(Xn,Xn,tau,K);
%DF = fdiff_fix_p_R0(Xn,tau,L,K);
DF = [DF(:,1),DF(:,4:end)];
F = function_F(Xn,Xn,tau,L,K);
while norm(F,1)>1E-14 & norm(F)<1E5
   Xn = Newton_step_fix_p_R0(Xn,F,DF);
   DF = function_DF(Xn,Xn,tau,K);
   DF = [DF(:,1),DF(:,4:end)];
   F = function_F(Xn,Xn,tau,L,K);
   norm(F,1)
end
% Clean data for real-valued solutions
Xn(1:7) = real(Xn(1:7));
Xn(10+4*N:10+4*N+N-1) = flip(conj(Xn(10+4*N+N+1:10+6*N)));
Xn(9+2*N:9+2*N+N-1) = flip(conj(Xn(9+2*N+N+1:9+4*N)));
Xn(8:8+N-1) = flip(conj(Xn(8+N+1:8+2*N)));
Xn(8+N) = real(Xn(8+N));
Xn(9+3*N) = real(Xn(9+3*N));
Xn(10+5*N) = real(Xn(10+5*N));
end