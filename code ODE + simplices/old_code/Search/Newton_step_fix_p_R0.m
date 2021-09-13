function Xn = Newton_step_fix_p_R0(X,F,DF)
Xn = [X(1);X(4:end)] - DF\F;
Xn = [Xn(1);X(2);X(3);Xn(2:end)];
end