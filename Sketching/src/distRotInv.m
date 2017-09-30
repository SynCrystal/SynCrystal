function d = distRotInv(v1,v2)
sz = size(v2);
vec = 1:sz(2);
mat = zeros(sz(1)*sz(2),sz(2));
for cnt = 1:sz(2)
    pos = circshift(vec,cnt-1,2); 
    temp = v2(:,pos)';
    mat(:,cnt) = temp(:);
end
v1 = v1';
v1 = v1(:);
v1 = repmat(v1,1,sz(2));
d = min(sqrt(sum((v1-mat).^2)))/sqrt(sz(1)*sz(2));