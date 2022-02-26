function norm= InfNuNorm_mat(mat)
% function norm= InfNuNorm_mat(mat)
% nu global variable

global nu
global use_intlab

if use_intlab && ~isintval(mat)
    mat = intval(mat);
end

m=(size(mat)-[1,1])/2;

m_row=m(1);m_col=m(2);

if m_row~=floor(m_row) || m_col~=floor(m_col)
    error('sizes of the matrix must be odd')
end

if ~use_intlab && ~isintval(mat)
    
    norm= max((nu.^abs(-m_row:m_row)*abs(mat)).*nu.^(-abs(-m_col:m_col)));%max((nu.^(-m_row:m_row)*abs(mat)).*nu.^(-(-m_col:m_col)));
else
    norm= max(sup((nu.^abs(-m_row:m_row)*abs(mat)).*nu.^(-abs(-m_col:m_col))));%max(sup((nu.^(-m_row:m_row)*abs(mat)).*nu.^(-(-m_col:m_col))));
end
% norm=-Inf;
% for k=-m_col:m_col
%     partial_norm=0;
%     for j=-m_row:m_row
%         if ~use_intlab
%             partial_norm=partial_norm+abs(mat(m_row+1+j,m_col+1+k))*nu^abs(j);
%         else
%             partial_norm=partial_norm+sup(abs(mat(m_row+1+j,m_col+1+k))*nu^abs(intval(j)));
%             
%         end
%     end
%     if ~use_intlab
%         partial_norm=partial_norm/(nu^abs(k));
%     else
%         partial_norm=partial_norm/sup(nu^abs(k));
%     end
%     norm=max(norm,partial_norm);
% end