function DS = derivative_saddle_fullop(xvw_PDE, G, Es, EsDelta)
% function DS = derivative_saddle_fullop(x_PDE, G, Es, EsDelta)
%
% INPUT
% x_PDE     PDE_vector with
%               size_real = 9
%               size_vector = 3
% G, Es     adequate length vectors for the linear equations
% EsDelta   same
% flagFullSize      DEFAULT 0, size of output operator
% OUTPUT
% DS        overoperator3

[x_PDE, v_PDE, w_PDE] = split(xvw_PDE);

DF = derivative_Hopf_fullop(x_PDE,G,Es);

D2Fv = second_der_Hopf_fullop(x_PDE,v_PDE);
D2Fw = second_der_Hopf_fullop(x_PDE,w_PDE);
D3Fvv = third_derHopf_fullop(x_PDE, v_PDE);


EsDeltavec = PDE_vector(EsDelta,x_PDE.size_real, x_PDE.size_vector, x_PDE.node_space);
EsDeltavec = reshape(EsDeltavec,x_PDE.node_space*2,x_PDE.node_time*2);
EsDelta = PDE_vec2vec(EsDeltavec);

x_PDE = reshape(x_PDE,x_PDE.node_space*2,x_PDE.node_time*2);

DsDF = zeros(length(x_PDE));
DsDF(1,:) = EsDelta(:).';
DsDF = compose_overoperator(DsDF, x_PDE);

Ds = overoperator3(0*DsDF, 0*DsDF, 0*DsDF,...
    DsDF, 0*DsDF, 0*DsDF,...
    0*DsDF, 2*DsDF , 0*DsDF);

Ds_fullop = full_operator(Ds);

DS = compose_fullop(DF, 0*DF, 0*DF,...
    D2Fv, DF, 0*DF,...
    D3Fvv+D2Fw, 2*D2Fv, DF) + Ds_fullop;