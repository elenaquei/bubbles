function [INTEGRAL,ERROR_MULTIPLE] = verifyquad_simpson(x,y,h)
% [INTEGRAL,ERROR_MULTIPLE] = verifyquad_simpson(x,y,h) 
% Computes an approximation (INTEGRAL) of the integral of a function f for
% which f(x0)=y0,...,f(xn)=yn, with x = [x0,...,xn] and y=[y0,...,yn] for
% column vectors yj, using Simpson's 1/3 rule, assuming xj = x0 + jh. Note
% that the vector of sample points x must have an odd number of terms,
% since the 1/3-rule uses an even number of subintervals. The error of this
% rule is ERROR_MULTIPLE*sup||f^(4)||.
n = size(x)*[0;1]-1;
d = size(y)*[1;0];
if even(n)==0
    error('input x must contain an odd number of evenly-spaced points.')
end
INTEGRAL = y(:,1) + y(:,n+1) + 2*sum(y(:,3:2:n-1),2) + 4*sum(y(:,2:2:n),2);
INTEGRAL = (h/3)*INTEGRAL;
ERROR_MULTIPLE = (x(end)-x(1))*h^4/180;
end