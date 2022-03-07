function zero_approx = Newton_handle(function_handle, approx, df)

max_iter = 100;
tol = 10^-12;
function_handle_t = @(t,x) function_handle(x);
for i = 1:max_iter
    if nargin<3
        DF = numjac(function_handle_t,0,approx,function_handle_t(0,approx),tol);
    else 
        DF = df(approx);
    end
    F = function_handle_t(0,approx);
    step = - DF \ F ;
    if i>1 && norm(step)<tol 
        return
    else
        %disp(norm(step))
    end
    zero_approx = approx + step;
    approx = zero_approx;
end
error('Newton did not converge')