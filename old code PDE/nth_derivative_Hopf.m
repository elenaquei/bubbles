function DnHopf = nth_derivative_Hopf(n, G, Es, x, varargin)
% function DnHopf = nth_derivative_Hopf(n, G, Es, x, varargin)
%
% INPUT
% n         integer, n<=4, order of derivative
% G         if n<=1, linear coefficients
% Es        if n<=1, linear coefficients of continuation equation
% x         PDE_vector
% varargin
%   v1...v(n-1)     directions of derivatives
%   flagFullSize    (DEFAULT 0) returns as big as possible
%   flagOp          (DEFAULT 1) returns an overoperator
% (if just one bool, interpreted as FullSize)
% OUTPUT
% DnHopf    overoperator or matrix, n-th derivative of the Hopf problem

if n<1 || n > 4
    error('Derivative impossible to compute')
end

if n==1 && isempty(G)
    error('Require linear infos for first derivative')
end

if ~isa(x,'PDE_vector')
    error('Required a PDE_vector as 4th input')
end

n_varagin = length(varargin);

if n_varagin<n-1
    error('For directional derivatives, directions required')
end
if n_varagin == n-1
    flagFullSize = 0;
    flagOp = 1;
elseif n_varagin == n
    flagFullSize = varargin{n_varagin};
    flagOp = 1;
else
    flagFullSize = varargin{n_varagin-1};
    flagOp = varargin{n_varagin};
end

switch n
    case 1
        DnHopf = derivative_Hopf(x, G, Es, flagFullSize, flagOp);
    case 2
        DnHopf = second_derHopf(x, varargin{1}, flagFullSize);
    case 3
        DnHopf = third_derHopf(x, varargin{1}, varargin{2},flagFullSize);
    case 4
        DnHopf = fourth_derHopf(x, varargin{1}, varargin{2}, ...
            varargin{3}, flagFullSize);
    otherwise
        error('Unrecognised derivative')
end
if flagOp && ~(n==1 && isempty(Es) )
	DnHopf = compose_overoperator (DnHopf, x);
end


