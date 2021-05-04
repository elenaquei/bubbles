function DH = derivative_saddle(xvw_PDE, G, Es, EsDelta, flagFullSize, flagOp)
% function DH = derivative_saddle(x_PDE, G, Es, EsDelta, flagFullSize, flagOp)
%
% INPUT
% x_PDE     PDE_vector with
%               size_real = 9
%               size_vector = 3
% G, Es     adequate length vectors for the linear equations
% EsDelta   same
% flagFullSize      DEFAULT 0, size of output operator
% flagOp            DEFAULT 0, output overoperator3
% OUTPUT
% DS        overoperator3

if nargin<5
    flagFullSize=0;
end

[x_PDE, v_PDE, w_PDE] = split(xvw_PDE);

DF = derivative_Hopf(x_PDE,G,Es,flagFullSize);

D2Fv = second_derHopf(x_PDE,v_PDE,flagFullSize);
D2Fw = second_derHopf(x_PDE,w_PDE,flagFullSize);
D3Fvv = third_derHopf(x_PDE, v_PDE, flagFullSize);

DsDF = 0*DF;
if flagFullSize
    EsDeltavec = PDE_vector(EsDelta,x_PDE.size_real, x_PDE.size_vector, x_PDE.node_space);
    EsDeltavec = reshape(EsDeltavec,x_PDE.node_space*2,x_PDE.node_time*2);
    EsDelta = PDE_vec2vec(EsDeltavec);
end
DsDF(1,:) = EsDelta(:).';

DH = [DF, 0*DF, 0*DF;
    DsDF+D2Fv, DF, 0*DF;
    D3Fvv+D2Fw, 2*DsDF+2*D2Fv, DF];
if nargin >5 && flagOp
    DH = compose_overoperator3(DH,xvw_PDE);
end