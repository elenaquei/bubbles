function [flag,Imin,Imax]=radii_polynomials(xBar,F_square,DFm, Am)
% function [flag,Imin,Imax]=radii_polynomials(xBar,F_square,DFm, Am)
%
% INPUT:
% xBar      Xi_vector, approximate solution of the problem;
% F_square  full_problem, vector field and scalar equations of the ODE
% DFm       complex matrix, derivative of F_square in xBar (DEFAULT:
%           derivative_to_matrix(derivative(F_square, xBar))
% Am        complex matrix, approximate inverse of DF (DEFAULT: inv(DF(xBar)) )
%
% OUTPUT:
% flag      0 if failed, 1 if successful
% Imin      positive real value, left bound of the interval;
% Imax      positive real value, right bound of the interval.
%
% In this function, the radii polynomial refering to the given problem and
% solution are computed.

global use_intlab
global Display

if nargin<3 || isempty(DFm)
    DFm=derivative_to_matrix(derivative(F_square,xBar,0));
end
if nargin<4 ||isempty(Am)
    Am=inv(DFm);
end
% change into intvals
if use_intlab
    xBar=intval(xBar);
    Am=intval(Am);
end

flag=0;
Imin=0;Imax=0;

Yvector=Y_bound(Am,xBar,F_square);

Z0vector=Z0_bound(DFm,Am,xBar);

Z1vector=Z1_bound(Am,xBar,F_square);

Z2vector=Z2_bound(Am,xBar,F_square);

p=@(r) Z0vector*r+Z1vector*r+Z2vector*r.^2-ones(size(Z0vector))*r+Yvector*ones(size(r));


if any(isnan(Z2vector+Z1vector+Z0vector+Yvector))
    return
end

[bool,Imin,Imax]=find_negative(Z2vector,Z1vector,Z0vector,Yvector);

if ~bool
    return
end

% plot
if Display
    xplot=0:0.000001:Imax*1.5;
    figure;
    pplot=p(xplot);
    plot(xplot,pplot,xplot,0*xplot,'k');
    hold on
    plot(Imin,0,'k*',Imax,0,'k*')
end
flag = 1;



return