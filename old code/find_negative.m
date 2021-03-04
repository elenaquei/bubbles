function [Imin,Imax]=find_negative(Z2vector,Z1vector,Z0vector,Yvector,Rstar)
% function [Imin,Imax]=find_negative(Z2vector,Z1vector,Z0vector,Yvector,Rstar)
% 
% this function computes the interval in which the radii polynomial is
% negative
%       P(r) = Yvector + (Z0vector + Z1vector -1) r + Z2vector r^2
%
% INPUT
% Z2vector              real vector, second order term
% Z1vector, Z0vector    real vector, first order term
% Yvector               real vector, zeroth order term
%                 all vectors must be of the same length
% Rstar                maximum radius
% OUTPUT
%  Imin, Imax           lower and upper bound of the interval in which the
%                       radii polynomial



if isa(Z2vector,'intval')
    a= sup(Z2vector);
    b = sup(Z1vector)+sup(Z0vector);
    b = b - 1;
    c = sup(Yvector);
    if any(sup(Z1vector+Z0vector)>1)
        error('Sum of Z1 and Z0 too big')
    end
else
    a = Z2vector;
    b = Z1vector + Z0vector -1;
    c = Yvector;
    if any((Z1vector+Z0vector)>1)
        error('Sum of Z1 and Z0 too big')
    end
end

Delta = b.^2 - 4*a.*c;

if any((Delta)<0)
    error('No interval found');
end

Imin = zeros(length(a),1);
Imax = zeros(length(a),1);

for i = 1:length(a)
    if a(i)>0
        Imin(i) = ( - b(i) - sqrt(Delta(i)))./(2*a(i));
        Imax(i) = ( - b(i) + sqrt(Delta(i)))./(2*a(i));
    else
        Imin(i) = -c(i)/b(i);
        Imax(i) = 10;
    end
end

Imin = max(Imin,0);
Imin = max(Imin);
Imax = min(Imax);
if Imin >Imax
    error('No interval found')
end

if Imin > Rstar
    error('Bounds bigger than Rstar')
else
    Imax = min(Imax, Rstar);
end

return
end