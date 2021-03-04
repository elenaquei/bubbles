function DDH = second_der_saddle(xvw_PDE, abc_PDE, flagFullSize)
% INPUT
% Xs    PDE_vector
% V     PDE_vector
% OUTPUT
% DDDHxDxD  overoperator

if nargin<3 || isempty(flagFullSize)
    flagFullSize = 0;
end

[x_PDE, v_PDE, w_PDE] = split(xvw_PDE);
[a_PDE, b_PDE, c_PDE] = split(abc_PDE);

D2Fa = second_derHopf(x_PDE,a_PDE,flagFullSize);
D2Fb = second_derHopf(x_PDE,b_PDE,flagFullSize);
D2Fc = second_derHopf(x_PDE,c_PDE,flagFullSize);

D3Fva = third_derHopf(x_PDE, v_PDE, a_PDE, flagFullSize);
D3Fvb = third_derHopf(x_PDE, v_PDE, b_PDE, flagFullSize);
D3Fwa = third_derHopf(x_PDE, w_PDE, a_PDE, flagFullSize);

D4Fvva = fourth_derHopf(x_PDE, v_PDE, v_PDE, a_PDE, flagFullSize);


DDH = overoperator3(D2Fa, 0*D2Fa, 0*D2Fa,...
    D3Fva+D2Fb, D2Fa, 0*D2Fa,...
    D4Fvva+D3Fwa+2*D3Fvb+D2Fc, 2*D3Fva+2*D2Fb, D2Fa);