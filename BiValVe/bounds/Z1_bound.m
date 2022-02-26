function [Z1vector,Z_mat,nodes_col,nodes_row]=Z1_bound(A,xBar,alpha)
% function Z1vector=Z1_bound(A,xBar,alpha)
%
% INPUT
% A      matrix, inverse of DF(x)
% xBar   Xi_vec, numerical solution
% alpha  full_problem
% OUTPUT
% Z1vector
% Z_mat
% nodes_col
% nodes_row


global use_intlab
global talkative
global Display

if ~isa(alpha,'full_problem')
    error('Wrong input')
end


N=xBar.size_vector; % N dimension of vector field 
M=xBar.size_scalar; % M number of scalar equations
F=2*xBar.nodes+1;   % F= (2n+1) number of nodes 
Deg=max(alpha.vector_field.deg_vector,2); % Deg degree of the vector field with a minimum of 2
nodes=xBar.nodes;   

DF =  derivative(alpha,xBar,1);    
DF.derivative_Fx_diagonal=0*DF.derivative_Fx_diagonal;
DFm = reshape(DF, nodes*(Deg-1)+1);
DHm = derivative_to_matrix(DFm);


nodes_new=DFm.nodes;
F_D_new=2*nodes_new+1;
diff=nodes_new-nodes;

screen=ones(size(DHm)); % this function will screen out the unwanted elements from DHm
screen(1:M,:)=0;

for i=0:N-1
    screen(M+i*F_D_new+diff+(1:F),...
        1:M)=0;
end

for i=0:N-1
    for j=0:N-1
        screen(M+i*F_D_new+diff+(1:F),...
            M+j*F_D_new+diff+(1:F))=0;
    end
end
if use_intlab
    screen=intval(screen);
end

DHm=screen.*DHm;

% at this point, the difference in size between DHm and A should be of 
% N*(Deg-2)*nodes. in order to multiply them correctly, we have to change
% the dimensions of A (remembering that on the diagonal we have the
% elements 1/(k alpha1 2pi)

% debugging check
if size(DHm,1)~=size(DHm,2)
    error('something got wrong :S');
elseif size(DHm,1) ~=  M+ N*F_D_new
    warning('probably something got wrong :-/');
end

A_big=extend_approximate_inverse(A,M,N, xBar.nodes,nodes_new);

% of DHm just the central columns need to be taken, getting to a matrix of
% the size 2(Deg-1)2nodes+3 X (Deg-1)2nodes+3
DHm_new=zeros(size(DHm,1), M+N*((Deg-1)*2*nodes+3));

if use_intlab
    DHm_new=intval(DHm_new);
end

DHm_new(:,1:M)=DHm(:,1:M);
col_new=(Deg-1)*2*nodes+3;
col_old=(size(DHm,2)-M)/N;
diff=(col_old-col_new)/2;
old_nodes = (col_old-1)/2;
new_nodes = (col_new -1)/2;

for col=1:N
    DHm_new(:,M+new_nodes+1+ col_new*(col-1)+(-old_nodes:old_nodes))=...
        DHm(:,M+old_nodes+1+ col_old*(col-1)+(-old_nodes:old_nodes));
end

Z_mat=A_big*DHm_new;

nodes_col=(Deg-1)*nodes+1;
nodes_row=((size(Z_mat,1)-M)/N-1)/2;
Z1vector=norm_Ximat(Z_mat,M,N,nodes_col,nodes_row);

% error handeling 
if any(Z1vector>1) %&& talkative>0
    fprintf('Z1 computed, %d\n',max(Z1vector));
    %error( 'Z1 is bigger than 1, no interval found' );
elseif talkative>2
   fprintf('Z1 computed, %d\n',max(Z1vector));
   fprintf('\n');
elseif talkative>1
   fprintf('Z1 computed, %d\n',max(Z1vector));
   fprintf('\n');
end

