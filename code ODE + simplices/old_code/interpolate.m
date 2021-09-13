function [X_s,PHI_s]= interpolate(DATA)
% DATA.x_i = [u_i;p_i;R_0_i;phi_i;epsilon_i;eta1;eta2;a0_i;a1_i;a2_i], i=0,1,2
% DATA.phi_i = [phi_1^i,phi_2^i], i=0,1,2
X_s = @(s) DATA.x_0 + s(1)*(DATA.x_1-DATA.x_0) + s(2)*(DATA.x_2-DATA.x_0);
PHI_s = @(s) DATA.phi_0 + s(1)*(DATA.phi_1-DATA.phi_0) + s(2)*(DATA.phi_2-DATA.phi_0);
end