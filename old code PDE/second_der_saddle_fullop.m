function DDH = second_der_saddle_fullop(xvw_PDE, abc_PDE)
% INPUT
% Xs    PDE_vector
% V     PDE_vector
% OUTPUT
% DDDHxDxD  overoperator

[x_PDE, v_PDE, w_PDE] = split(xvw_PDE);
[a_PDE, b_PDE, c_PDE] = split(abc_PDE);

D2Fa = second_der_Hopf_fullop(x_PDE,a_PDE);
D2Fb = second_der_Hopf_fullop(x_PDE,b_PDE);
D2Fc = second_der_Hopf_fullop(x_PDE,c_PDE);

D3Fva = third_derHopf_fullop(x_PDE, v_PDE, a_PDE);
D3Fvb = third_derHopf_fullop(x_PDE, v_PDE, b_PDE);
D3Fwa = third_derHopf_fullop(x_PDE, w_PDE, a_PDE);

D4Fvva = fourth_derHopf_fullop(x_PDE, v_PDE, v_PDE, a_PDE);


DDH = compose_fullop(D2Fa, 0*D2Fa, 0*D2Fa,...
    D3Fva+D2Fb, D2Fa, 0*D2Fa,...
    D4Fvva+D3Fwa+2*D3Fvb+D2Fc, 2*D3Fva+2*D2Fb, D2Fa);