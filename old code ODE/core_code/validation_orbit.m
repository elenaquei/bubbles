function [bool,Imin,Imax] = validation_orbit(Func, xXi_vec,DF_mat, A)
% function bool = validation_orbit(Func, xXi_vec)
% 
% validates the pointwise orbit proposed with respect to the given problem
% Automatically switches on intlab
%
% INPUT
% Func     instance of full_problem
% xXi_vec  instance of Xi_vector, proposed numerical solution
% OUTPUT
% bool     1 if validation successful, 0 otherwise
global use_intlab
global talkative
global RAD_MAX

temp_intlab = use_intlab;
if isempty(talkative)
    talkative = 0;
end

bool = 0;
Imin = [];
Imax = [];

if nargin<3
    DF =  derivative(Func,xXi_vec,0);
    DF_mat = derivative_to_matrix(DF);
end
if nargin<4
    A  = inv(DF_mat);
end

use_intlab = 1;

Y_vector = Y_bound_new(A,xXi_vec,Func);
Z0_vector=Z0_bound(DF_mat,A,xXi_vec);
if any(Z0_vector>1)
    return
end
Z1_vector=Z1_bound_new(A,xXi_vec,Func);
if any(Z1_vector+Z0_vector>1)
    return
end
Z2_vector= Z2_bound_new(A,xXi_vec,Func);
[bool,Imin,Imax]=find_negative(Z2_vector,Z1_vector,Z0_vector,Y_vector);

if ~bool
    return
end
% check
if Imax>RAD_MAX
    Imax = RAD_MAX;
end
if Imin >RAD_MAX
    bool = 0;
    Imin = [];
    Imax = [];
    display('RAD_MAX too small')
    return
end

if talkative
    fprintf('\n The interval is [ %e, %e ].\n\n',Imin,Imax);
end

use_intlab = temp_intlab;
bool = 1;

end