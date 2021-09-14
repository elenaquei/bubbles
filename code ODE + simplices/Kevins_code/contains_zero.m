function yn = contains_zero(u)
% Check if real INTVAL u contains zero
if ~isintval(u)
    u=intval(u);
end
yn = (inf(u)<=0)*(0<=sup(u));
end

