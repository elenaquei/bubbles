function G_Hopf_new = new_amplitude(X1)
%function G_Hopf_new = new_amplitude(G_Hopf, X0, X1)
% 
% function setting G2 to <K^2 hat z, z > = 1
% where hat z is the average between z0 and z1


y = X1.vector1D;
z = squeeze(X1.vector2D);
a = X1.parameters(end);

K = -X1.node_space:X1.node_space;
J = -X1.node_time:X1.node_time;

K2 = repmat(K,2*X1.node_time+1,1);
J2 = repmat(J(:), 1, 2*X1.node_space+1,1);

i_J_z = 1i*J2.*z;
G1 = conj([0,0,0,0*y(:).',i_J_z(:).']);

K_y = 1i*y.*K;
G3 = conj([0,0,0,K_y(:).',0*z(:).']);

%G2 = [0,0,0*y_hat(:).',conj(z_hat(:)).']; %old
J2z = (J2.^2.*z);
G2 = conj([0,0,0,0*y(:).',(J2z(:)).']); % NEW 

y_2D = embedding(y,z);
y_plus_az = y_2D+a*z;
G4 = conj([0,0,0, 0*y(:).', 1i*K2(:).'.*y_plus_az(:).']);

G_Hopf_new= [G1;G2;G3;G4];
% const_G_Hopf = [const_G1,const_G2,const_G3,const_G4];


return
a_hat = (X0.parameters(3)*2+0*X1.parameters(3))/2;
y_hat = (X0.vector1D*2+0*X1.vector1D)/2;
z_hat = squeeze((X0.vector2D *2+0* X1.vector2D)/2);

K = -X0.node_space:X0.node_space;
J = -X0.node_time:X0.node_time;

K2 = repmat(K,2*X0.node_time+1,1);
J2 = repmat(J(:), 1, 2*X0.node_space+1,1);

i_J_z_hat = 1i*J2.*z_hat;
G1 = conj([0,0,0,0*y_hat(:).',i_J_z_hat(:).']);
%const_G1 =0;
K_y_hat = 1i*y_hat(:).*K(:);
G3 = conj([0,0,0,K_y_hat(:).',0*z_hat(:).']);
%const_G3 = 0;
%G2 = [0,0,0*y_hat(:).',conj(z_hat(:)).']; %old
J2z = (J2.^2.*z_hat);
G2 = conj([0,0,0,0*y_hat(:).',(J2z(:)).']); % NEW 
%const_G2 = -1;
y_2D = embedding(y_hat,z_hat);
y_plus_az = y_2D+a_hat*z_hat;
G4 = conj([0,0,0, 0*y_hat(:).', 1i*K2(:).'.*y_plus_az(:).']);
%const_G4 = 0;

G_Hopf_new= [G1;G2;G3;G4];
%const_G_Hopf = [const_G1,const_G2,const_G3,const_G4];
% constants do not change over the iterations

end