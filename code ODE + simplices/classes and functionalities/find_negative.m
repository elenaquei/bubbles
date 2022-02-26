function [bool,Imin,Imax]=find_negative(Z2vector,Z1vector,Z0vector,Yvector)
% function [bol,Imin,Imax]=find_negative(Z2vector,Z1vector,Z0vector,Yvector)
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
% OUTPUT
% bool                  success boolean
%  Imin, Imax           lower and upper bound of the interval in which the
%                       radii polynomial
global use_intlab
bool = 0;
Imin = [];
Imax = [];

if use_intlab
    a= intval(Z2vector)+1;
    b = intval(Z1vector)+intval(Z0vector);
    b = b - 1;
    c = intval(Yvector);
    if any(sup(Z1vector+Z0vector)>1)
        warning('Sum of Z1 and Z0 too big')
        return
    end
else
    use_intval = 0;
    a = Z2vector;
    b = Z1vector + Z0vector -1;
    c = Yvector;
    if any((Z1vector+Z0vector)>1)
        warning('Sum of Z1 and Z0 too big')
        return
    end
end

Delta = b.^2 - 4*a.*c;

if ( use_intlab == 1 && any(inf(Delta)<0) )||( use_intlab == 0 && any((Delta)<0) )
    error('VAL:NegDelta','No interval found');
end

Imin = zeros(length(a),1);
Imax = zeros(length(a),1);
if use_intlab 
    Imin = intval(Imin);
    Imax = intval(Imax);
end

for i = 1:length(a)
    if a(i)>0
        Imin(i) = ( - b(i) - sqrt(Delta(i)))./(2*a(i));
        Imax(i) = ( - b(i) + sqrt(Delta(i)))./(2*a(i));
    else
        Imin(i) = -c(i)/b(i);
        Imax(i) = 10;
    end
end
if use_intlab 
    Imin = sup(Imin);
    Imax = inf(Imax);
end
Imin = max(Imin,0);

Imin = max(Imin);
Imax = min(Imax);
if Imin >Imax
    Imin = Imax;
end

bool = 1;

return
end
