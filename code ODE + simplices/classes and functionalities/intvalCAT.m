function c=intvalCAT(dim,a,b)
% cat(DIM,A,B) concatenates the intval arrays A and B along
%     the dimension DIM.  
%     cat(2,A,B) is the same as [A,B].
%     cat(1,A,B) is the same as [A;B].
%  
%     B = cat(DIM,A1,A2,A3,A4,...) concatenates the input
%     arrays A1, A2, etc. along the dimension DIM.

if ~isintval(a) || ~isintval(b)
    warning('this function is just for intvals')
end

if size(dim,2)~=1 || size(dim,1)~=1  || any(size(size(dim))~= [1,2])
    error('firt input must be an integer')
end

size_A=size(a);
size_B=size(b);

if length(size_A)~=length(size_B)
    error('second and third inputs must have the same number of dimensions');
end

if any(size_A(1:dim-1)~= size_B(1:dim-1)) || any(size_B(dim+1:end)~= size_B(dim+1:end))
    error('sizes must agree')
end

index=zeros(1,length(size_A));index(dim)=1;

size_C=size_A+index.*size_B;
c=intval(zeros(size_C));
switch length(size_A)
    case 2
        if dim==1
            c(1:size_A(1),1:size_A(2))=a(:,:);
            c(size_A(1)+1:end,1:size_B(2))=b(:,:);
        elseif dim==2
            c(1:size_A(1),1:size_A(2))=a(:,:);
            c(1:size_A(1),size_A(2)+1:end)=b(:,:);
        else
            error('dim is not consistent');
        end
    case 3
        if dim==1
            c(1:size_A(1),1:size_A(2),1:size_A(3))=a(:,:,:);
            c(size_A(1)+1:end,1:size_B(2),1:size_C(3))=b(:,:,:);
        elseif dim==2
            c(1:size_A(1),1:size_A(2),1:size_A(3))=a(:,:,:);
            c(1:size_B(1),size_A(2)+1:end,1:size_B(3))=b(:,:,:);
        elseif dim==3
            c(1:size_A(1),1:size_A(2),1:size_A(3))=a(:,:,:);
            c(1:size_B(1),1:size_B(2),size_A(3)+1:end)=b(:,:,:);
        else
            error('dim is not consistent');
        end
    case 4
        if dim==1
            c(1:size_A(1),1:size_A(2),1:size_A(3),1:size_A(4))=a(:,:,:,:);
            c(size_A(1)+1:end,1:size_B(2),1:size_B(3),1:size_B(4))=b(:,:,:,:);
        elseif dim==2
            c(1:size_A(1),1:size_A(2),1:size_A(3),1:size_A(4))=a(:,:,:,:);
            c(1:size_B(1),size_A(2)+1:end,1:size_B(3),1:size_B(4))=b(:,:,:,:);
        elseif dim==3
            c(1:size_A(1),1:size_A(2),1:size_A(3),1:size_A(4))=a(:,:,:,:);
            c(1:size_B(1),1:size_B(2),size_A(3)+1:end,1:size_B(4))=b(:,:,:,:);
        elseif dim==4
            c(1:size_A(1),1:size_A(2),1:size_A(3),1:size_A(4))=a(:,:,:,:);
            c(1:size_B(1),1:size_B(2),1:size_B(3):end,size_A(4)+1:end)=b(:,:,:,:);
        else
            error('dim is not consistent');
        end
        
    case 5
        if dim==1
            c(1:size_A(1),1:size_A(2),1:size_A(3),1:size_A(4),1:size_A(5))=a(:,:,:,:,:);
            c(size_A(1)+1:end,1:size_B(2),1:size_B(3),1:size_B(4),1:size_B(5))=b(:,:,:,:,:);
        elseif dim==2
            c(1:size_A(1),1:size_A(2),1:size_A(3),1:size_A(4),1:size_A(5))=a(:,:,:,:,:);
            c(1:size_B(1),size_A(2)+1:end,1:size_B(3),1:size_B(4),1:size_B(5))=b(:,:,:,:,:);
        elseif dim==3
            c(1:size_A(1),1:size_A(2),1:size_A(3),1:size_A(4),1:size_A(5))=a(:,:,:,:,:);
            c(1:size_B(1),1:size_B(2),size_A(3)+1:end,1:size_B(4),1:size_B(5))=b(:,:,:,:,:);
        elseif dim==4
            c(1:size_A(1),1:size_A(2),1:size_A(3),1:size_A(4),1:size_A(5))=a(:,:,:,:,:);
            c(1:size_B(1),1:size_B(2),1:size_B(3):end,size_A(4)+1:end,1:size_A(5))=b(:,:,:,:,:);
        elseif dim==5
            c(1:size_A(1),1:size_A(2),1:size_A(3),1:size_A(4),1:size_A(5))=a(:,:,:,:,:);
            c(1:size_B(1),1:size_B(2),1:size_B(3):end,1:size_B(4),size_A(5)+1:end)=b(:,:,:,:,:);
        else
            error('dim is not consistent');
        end
        
    otherwise
        error('such a high dimension is not possible at the moment')
        
end


% function c = vertcat(varargin)
% %VERTCAT      Implements  [a(1) ; a(2) ; ...]  for Taylor
% %
% 
% % written  05/21/09     S.M. Rump
% %
% 
%   a = taylor(varargin{1});
%   c.size = a.size;
%   index = reshape(1:prod(a.size),a.size)';
%   c.t = a.t(:,index);             % transposed of first element a
%   cols = c.size(2);
% 
%   for i=2:length(varargin)
%     a = taylor(varargin{i});
%     if cols~=a.size(2)
%       error('dimension do not fit')
%     end
%     c.size(1) = c.size(1) + a.size(1);
%     index = reshape(1:prod(a.size),a.size)';
%     c.t = [ c.t a.t(:,index) ];   % horzcat for columnwise stored arrays
%   end
%   
%   index = reshape(1:prod(c.size),fliplr(c.size))';
%   c.t = c.t(:,index);             % transpose result
% 
%   c = class(c,'taylor');