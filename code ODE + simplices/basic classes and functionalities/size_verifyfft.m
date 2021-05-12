function size_z = size_verifyfft(z)
% function size_z = size_verifyfft(z)
% 
% returns the size of verifyfft_in(z)
len=length(z);
next_pow=ceil(log2(len));
size_z = 2^next_pow;