function conv_u_v = conv2INTVAL(u,v, option)
% function conv2INTVAL(u,v)
% 
% u,v       2D vectors intval

if ~isa(v,'intval')
    v = intval(v);
end
if ~isa(u,'intval')
    u = intval(u);
end

size_u = size(u);
size_v = size(v);

max_size = max(2*size_u-1,2*size_v-1);

fft_size = 2 .^(ceil(log2(max_size)));

u_big = intval(zeros(fft_size));
v_big = intval(zeros(fft_size));

u_big(1:size_u(1),1:size_u(2)) = u;
v_big(1:size_v(1),1:size_v(2)) = v;

conv_u_v = altifftn(altfftn(u_big).*altfftn(v_big));

if all(option == 'same')
    max_size = max(size_u,size_v);
end

conv_u_v = conv_u_v(1:max_size(1),1:max_size(2));

% warning('To update')
% conv_u_v = real_conv(real(u),real(v),option) + ...
%     1i*( real_conv(real(u),imag(v),option) + real_conv(imag(u),real(v),option) )...
%     - real_conv(imag(u),imag(v),option);
% 

end

function conv_u_v_real = real_conv(u,v,option)
%warning('using conv2 is not verified \n(and has non-deterministic behavior in Matlab2014b)\nUSE VERIFYFFT')
conv_u_v_infinf = conv2(inf(u),inf(v),option);
conv_u_v_supsup = conv2(sup(u),sup(v),option);
conv_u_v_infsup = conv2(inf(u),sup(v),option);
conv_u_v_supinf = conv2(sup(u),inf(v),option);

conv_u_v_sup = max(max(conv_u_v_infinf, conv_u_v_infsup), max(conv_u_v_supsup, conv_u_v_supinf));
conv_u_v_inf = min(min(conv_u_v_infinf, conv_u_v_infsup), min(conv_u_v_supsup, conv_u_v_supinf));

conv_u_v_real = infsup(conv_u_v_inf,conv_u_v_sup);
end
