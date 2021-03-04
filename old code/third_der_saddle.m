function DDH = third_der_saddle(xvw_PDE, abc_PDE, flagFullSize)
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

D3Faa = third_derHopf(x_PDE, a_PDE, a_PDE, flagFullSize);
D3Fab = third_derHopf(x_PDE, a_PDE, b_PDE, flagFullSize);
D3Fac = third_derHopf(x_PDE, a_PDE, c_PDE, flagFullSize);
D3Fbb = third_derHopf(x_PDE, b_PDE, b_PDE, flagFullSize);

D4Fvaa = fourth_derHopf(x_PDE, a_PDE, a_PDE, v_PDE, flagFullSize);
D4Fwaa = fourth_derHopf(x_PDE, a_PDE, a_PDE, w_PDE, flagFullSize);
D4Fvab = fourth_derHopf(x_PDE, v_PDE, b_PDE, a_PDE, flagFullSize);


DDH = overoperator3(D3Faa, 0*D3Faa, 0*D3Faa,...
    D4Fvaa+D3Fab, D3Faa, 0*D3Faa,...
    D4Fwaa+2*D4Fvab+D3Fbb+2*D3Fac, 2*(D4Fvaa+D3Fab), D3Faa);