function C = cswft2_fwd(x,h,sigma)
% cswft2_fwd compute the compactly supported windowed Fourier transform using
% Gaussian window
%
% Inputs
%   x           sz1-by-sz2 matrix
%
% Optional Inputs
%   h          h(1) is the spacing of samples in the vertical direction in the image
%              h(2) is the spacing of samples in the horizontal direction in the image
%   sigma       width parameter of the Gaussian window
%
% Outputs
%   C           Cell array of windowed Fourier coefficients.
%               C(k1,k2,x1,x2) is the coefficient at
%                   - position (x1,x2)
%                   - frequency (k1,k2)
%
%by Haizhao Yang

[sz1,sz2] = size(x);
if nargin < 3, sigma =max(10,min(sz1,sz2)/10); end;
if nargin < 2, h = [1, 1]; end;

splp1 = ceil(sz1/h(1));
splp2 = ceil(sz2/h(2));
C = zeros(6*sigma,6*sigma,splp1,splp2);
for cnt1 = 1:splp1
    for cnt2 = 1:splp2
        ct1 = (cnt1-1)*h(1)+1; ct2 = (cnt2-1)*h(2)+1;
        idx1 = ct1-sigma*3:ct1+sigma*3-1;
        pos1 = mod(idx1,sz1)+1;
        temp = x(pos1,:);
        idx2 = ct2-sigma*3:ct2+sigma*3-1;
        pos2 = mod(idx2,sz2)+1;
        temp = temp(:,pos2);
        temp = temp.*gau2d(idx1,idx2,ct1,ct2,sigma);
        C(:,:,cnt1,cnt2) = fftshift(fft2(ifftshift(temp)))/sigma/6;
    end
end
end

function g = gau2d(x,y,xc,yc,sigma)
g = (((x'-xc).^2)*ones(1,length(y)) + ones(length(x),1)*((y-yc).^2))./(2*sigma^2);
g = exp(-g);
end