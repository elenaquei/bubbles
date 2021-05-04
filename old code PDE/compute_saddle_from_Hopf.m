function xvw = compute_saddle_from_Hopf(x_PDE,G,Es0,Es1,contEs0, contEs1)
% function xvw = compute_saddle_from_Hopf(x_PDE,G,Es0,Es1,contEs0, contEs1)
if x_PDE.size_real~=3 || x_PDE.size_vector~=1
    error('Input size unexpected')
end

DxF = derivative_Hopf(x_PDE,G,Es0);
DsF = 0*DxF(:,1);

C_Es = compute_linear_func(PDE_vec2vec(x_PDE), Es1-Es0, contEs1 - contEs0);
DsF(1) = C_Es;

x_prime = -DxF\DsF ;

% test = DxF*x_prime.' + DsF;

DxxF_xprime_xprime = second_derHopf(x_PDE,PDE_vector(x_prime,x_PDE))*x_prime;
DsX = 0*DsF;
DsX(1) =  compute_linear_func(x_prime.',Es1-Es0,0);

x_primeprime = - DxF\(DxxF_xprime_xprime + 2*DsX);

xvw = merge(x_PDE, PDE_vector(x_prime,x_PDE), PDE_vector(x_primeprime,x_PDE));