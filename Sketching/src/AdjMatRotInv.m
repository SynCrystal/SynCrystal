function mat = AdjMatRotInv(dataIn,sigma)
sz = size(dataIn);
mat = zeros(sz(1),sz(1));
for cnt1 = 1:sz(1)-1
    for cnt2 = cnt1+1:sz(1)
        v1 = reshape(dataIn(cnt1,:,:),[sz(2),sz(3)]);
        v2 = reshape(dataIn(cnt2,:,:),[sz(2),sz(3)]);
        mat(cnt1,cnt2) = gauss(distRotInv(v1,v2),sigma);
    end
end
mat = mat'+ mat + eye(sz(1));

