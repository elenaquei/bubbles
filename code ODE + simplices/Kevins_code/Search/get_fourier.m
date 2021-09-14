function a = get_fourier(s,f,phi,N,h)
a = zeros(2*N+1,1);
for k=-N:N
    a(k+N+1) = (phi/2/pi)*verifyquad_simpson(s,exp(-1i*phi*k*s).*f,h);
end
end