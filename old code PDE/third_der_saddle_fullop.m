function DDH = third_der_saddle_fullop(xvw_PDE, abc_PDE, ~)
% INPUT
% Xs    PDE_vector
% V     PDE_vector
% OUTPUT
% DDDHxDxD  full_operator (3)

[x_PDE, v_PDE, w_PDE] = split(xvw_PDE);
[a_PDE, b_PDE, c_PDE] = split(abc_PDE);

D3Faa = third_derHopf_fullop(x_PDE, a_PDE, a_PDE);
D3Fab = third_derHopf_fullop(x_PDE, a_PDE, b_PDE);
D3Fac = third_derHopf_fullop(x_PDE, a_PDE, c_PDE);
D3Fbb = third_derHopf_fullop(x_PDE, b_PDE, b_PDE);

D4Fvaa = fourth_derHopf_fullop(x_PDE, a_PDE, a_PDE, v_PDE);
D4Fwaa = fourth_derHopf_fullop(x_PDE, a_PDE, a_PDE, w_PDE);
D4Fvab = fourth_derHopf_fullop(x_PDE, v_PDE, b_PDE, a_PDE);


DDH = compose_fullop(D3Faa, 0*D3Faa, 0*D3Faa,...
    D4Fvaa+D3Fab, D3Faa, 0*D3Faa,...
    D4Fwaa+2*D4Fvab+D3Fbb+2*D3Fac, 2*(D4Fvaa+D3Fab), D3Faa);