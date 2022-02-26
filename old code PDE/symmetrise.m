function x_sym = symmetrise(x)
x_sym = x;

x_sym.parameters = (x.parameters + conj(x.parameters))/2;

y_sym = (x.vector1D(:,:) + conj(x.vector1D(:,end:-1:1)))/2;

z_sym = (x.vector2D(:,:,:) + conj(x.vector2D(:,end:-1:1,end:-1:1)))/2;

x_sym.vector1D = y_sym;

x_sym.vector2D = z_sym;
end