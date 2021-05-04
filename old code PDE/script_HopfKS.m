global nu
global display
display = 0;
nu = [1.1,1.1];

% ordering of parameters: nu, lambda, a 

try 
    intval(1);
catch
    addpath(genpath('./'));
    startintlab
end

fprintf('Code started at %s\n',datestr(now,13))

n_iter = 4;
step_size = 10^-6;
Rstar = 10^-5; %upper bound of the radius
% single point Hopf validation works with
% node_space = 30;
% node_time = 20;

node_space = 30;
node_time = 9;

[X, G_Hopf, const_G_Hopf, Es_Hopf, const_Es_Hopf]= Hopf_first_step(node_space, node_time);

X = symmetrise(X);

A = matrix_AHopf(X,G_Hopf,Es_Hopf,1);
A_dagger =  derivative_Hopf(X,G_Hopf,Es_Hopf);

%validate starting point
% Y_bound=YboundHopf(A,X,G_Hopf,const_G_Hopf, Es_Hopf,const_Es_Hopf);
% 
% Z0 = Z0Hopf(A_dagger,X,G_Hopf,Es_Hopf);
% 
% Z1 = Z1Hopf(A,A_dagger,X,G_Hopf,Es_Hopf);
% 
% Z2 = Z2Hopf(A,X,Rstar);
% try
%    [Imin,Imax]=find_negative(Z2,Z1,Z0,Y_bound,Rstar);
% catch
%    error('starting point not validated')
% end
% %return 
% - % - % - % - % - % - % - % - % - % - % - % - % - % - %
X0 = X;
Es0 = Es_Hopf;
const_Es0 = const_Es_Hopf;
G0 = G_Hopf;
const_G0 = const_G_Hopf; % never updated !
A0 = A;
A_dagger0 = sparse(A_dagger);
clear A_dagger A

a0 = X0.parameters(3);
Imin_start = 10^-9;% Imin

bound_Y = zeros(5, n_iter);
bound_Z0 = zeros(5, n_iter);
bound_Z1 = zeros(5, n_iter);
bound_Z2 = zeros(5, n_iter);
Imin_vec = zeros(1, n_iter);
Imax_vec = zeros(1, n_iter);

Xnorm_vec = zeros(5, n_iter);

for iter = 1:n_iter
    
    [X1,Es1,const_Es1, G1] = archlength(X0, G0, const_G_Hopf, Es0, const_Es0, step_size);
    
    X1 = symmetrise(X1);
    
    A1 = matrix_AHopf(X1,G1,Es1,1);
    A_dagger1 =  sparse(derivative_Hopf(X1,G1,Es1));
    
    Z0cont = Z0HopfCont(A_dagger0,A_dagger1,X0,X1);
    fprintf('Z0cont computed, %e\n',max(Z0cont));
    % works
    
    if any(Z0cont>1)
        error('G2 still wrong')
    end
    
    Z1cont = Z1HopfCont(A0, A1, A_dagger0,A_dagger1,X0,X1);
    fprintf('Z1cont computed, %f\n',max(Z1cont));
    % works
    
    if any((Z1cont+Z0cont)>1)
        error('G2 still wrong')
    end
    
    Z2cont = Z2HopfCont(A0, A1, X0, X1, Rstar);
    fprintf('Z2cont computed, %e\n',max(Z2cont));
    % works - and much simpler
    
    
    Ycont = YHopfCont(A0, A1,X0,X1,G0,G1,const_G_Hopf,...
        Es0,Es1, const_Es0, const_Es1);
    fprintf('Ycont computed, %e\n',max(Ycont));
    % works % WAY TOO SLOW!
    % all due to derivative_Hopf(Xs, G, Es_s, 1)
    % check it out!
    
    bound_Y(:,iter) = Ycont;
    bound_Z0(:,iter) = Z0cont;
    bound_Z1(:,iter) = Z1cont;
    bound_Z2(:,iter) = Z2cont;
    
    [Imin,Imax]=find_negative(Z2cont,Z1cont,Z0cont,Ycont,Rstar);
    
    Imin_vec(iter) = Imin;
    Imax_vec(iter) = Imax;
    Xnorm_vec(:,iter) = norm(X0);
    
    % check if Hopf was in this segment
    a1 = X1.parameters(3);
    if a0*a1<0  
        disp('EUREKA! we went throught the Hopf bifurcation')
        if midrad(a0,Imin_start)*midrad(a1,Imin)<0
            disp('EUREKA! we validated the Hopf bifurcation')
            Hop_validated = 1;
            %break
        end
    end
    
    % update
    X0 = X1;
    Es0 = Es1;
    const_Es0 = const_Es1;
    A0 = A1;
    A_dagger0 = A_dagger1;
    G0 = G1;
    a0 = a1;
    fprintf('Iteration %i completed at %s\n',iter,datestr(now,13))
end
save('longer_valdation_Hopf_1000iter')
plot(X,'Linewidth',4)
plot(X1,'Linewidth',4)
plot(X-X1,'Linewidth',4)

return
% the saddle problem will be solved later

% saddle problem
X0_saddle = compute_saddle_from_Hopf(X0,G_Hopf,Es0_Hopf,Es1_Hopf,const_Es0, const_Es1); 
% works

EsDelta = Es1_Hopf-Es0_Hopf;
const_EsDelta = const_Es1-const_Es0;
y = evaluate_saddle(X0_saddle, G_Hopf, const_G_Hopf,Es_Hopf,const_Es_Hopf,EsDelta,const_EsDelta);
% works

[X1_saddle,G_saddle, const_G_saddle,Es1_saddle,const_Es1_saddle] = archlength_saddle(X0_saddle, ...
    G_Hopf, const_G_Hopf, Es_Hopf, const_Es_Hopf,EsDelta,const_EsDelta, step_size);
y = evaluate_saddle(X1_saddle, G_Hopf, const_G_Hopf,Es1_Hopf,const_Es1_saddle,EsDelta,const_EsDelta);

X0 = X0_saddle;
X1 = X1_saddle;
Es0 = [Es_Hopf,0*Es_Hopf,0*Es_Hopf];
const_Es0 = const_Es_Hopf;
A_dagger0 = derivative_saddle(X0,G_Hopf,Es_Hopf,EsDelta);
A0 = matrix_Asaddle(A_dagger0,X0,[],[],1); % including A_big_saddle

A_dagger1 = derivative_saddle(X1,G_Hopf,Es_Hopf,EsDelta);
A1 = matrix_Asaddle(A_dagger1,X1,[],[],1); % including A_big_saddle

% makes sense for A and A_dagger to have different sizes

Ycont = YsaddleCont(A0, A1,X0,X1,G_Hopf,const_G_Hopf,Es0_Hopf,Es1_Hopf, const_Es0, const_Es1,EsDelta, const_EsDelta);

Z0cont = Z0saddleCont(A_dagger0,A_dagger1,X0,X1,G,Es0_Hopf,Es1_Hopf);

Z1cont = Z1saddleCont(A0, A1, A_dagger0,A_dagger1,X0,X1,G,Es0_Hopf,Es1_Hopf);

Z2cont = Z2saddleCont(A0, A1, X0, X1);

