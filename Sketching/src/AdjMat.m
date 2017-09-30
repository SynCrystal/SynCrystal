function mat = AdjMatRotInv(dataIn,sigma)
sz = size(dataIn);
mat = zeros(sz(3),sz(3));
for cnt1 = 1:sz(3)-1
    for cnt2 = cnt1+1:sz(3)
        mat(cnt1,cnt2) = gauss(norm(dataIn(:,:,cnt1)-dataIn(:,:,cnt2))/sqrt(sz(1)*sz(2)),sigma);
    end
end
mat = mat'+ mat + eye(sz(3));

