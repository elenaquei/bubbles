function AB = infsup_help(A,B)
% function AB = infsup_help(A,B)
%
% accepts any two matrices as input and creates the corresponding interval
% matrix

AB = A;

matrix_inf_real = min(real(A ), real(B ));
matrix_sup_real = max(real(A ), real(B ));

matrix_inf_imag = min(imag(A ), imag(B ));
matrix_sup_imag = max(imag(A ), imag(B ));

matrix_inf = matrix_inf_real + 1i*matrix_inf_imag;
matrix_sup = matrix_sup_real + 1i*matrix_sup_imag;

AB = infsup(matrix_inf, matrix_sup);