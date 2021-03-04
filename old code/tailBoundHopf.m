function B = tailBoundHopf(n,nu, lambda, nmodes)
% function B = tailBoundHopf(n,nu, lambda, nmodes)
% 
% computes a bound for the tail of
%
% K_2^n/(lambda*(K_2^4 nu - K_2^2) )
%
% startig a nmodes 

if nmodes<=0
    error('This is a tail bound, number of modes needs to be positive')
end

f=@(x) x.^n./(lambda*(x.^4 * nu - x.^2) );

if nu >0
    argmax = sqrt(1/abs(nu));
else
    argmax = 0;
end

if n<4
    k_2 = (2-n)/(4*(nu*(4-n)));
    if k_2>0
        if isa(k_2,'intval')
            argmax = max(sup(argmax),sup(sqrt(k_2)));
        else
            argmax = max(argmax,sqrt(k_2));
        end
    end
end

if argmax < nmodes
    if isa(nu, 'intval')
        B = sup(abs(f(nmodes)));
    else
        B = abs(f(nmodes));
    end
else
    if isa(nu, 'intval')
        %vec_argmax = floor(argmax):ceil(argmax);
        B = max(sup(abs(f(nmodes:ceil(argmax)))));
    else
        B = max((abs(f(nmodes:ceil(argmax)))));%max(abs(f(nmodes)), abs(f(floor(argmax))), abs(f(ceil(argmax))));
    end
end

if n == 4 % if same power, also consider limit to infinity 
    if isa(nu, 'intval')
        B = max( B, sup(abs(1/(nu*lambda))));
    else
        B = max( B, abs(1/(nu*lambda)));
    end
end