function Z2cont = Z2saddleCont(A0_in, A1_in, X0, X1)

A0 = fullopA_saddle(A0_in,X0);
A1 = fullopA_saddle(A1_in,X1);

Xs = intval(X0,X1);
XDelta = X0-X1;

As = intval(A0,A1);

Z2cont = norm(As*second_der_saddle_fullop(Xs, XDelta));