function C = cswft1_fwd(x,h,sigma)
% cswft1_fwd compute the 1D compactly supported windowed Fourier transform using
% Gaussian window
%
% Inputs
%   x          a vector of size N
%
% Optional Inputs
%   h          h is the spacing of samples
%   sigma       width parameter of the Gaussian window
%
% Outputs
%   C           Cell array of windowed Fourier coefficients.
%               C(k,x) is the coefficient at
%                   - position x
%                   - frequency k
%
%by Haizhao Yang

N = length(x);
if nargin < 3, sigma =max(10,N/10); end;
if nargin < 2, h = 1; end;

splp1 = ceil(N/h);
C = zeros(6*sigma,splp1);
for cnt1 = 1:splp1
    ct1 = (cnt1-1)*h+1;
    idx1 = ct1-sigma*3:ct1+sigma*3-1;
    pos1 = mod(idx1,N)+1;
    temp = x(pos1);
    temp = temp.*gau1d(idx1,ct1,sigma);
    C(:,cnt1) = fftshift(fft(ifftshift(temp)))/sqrt(sigma/6);
end
end

function g = gau1d(x,xc,sigma)
g = (((x-xc).^2))./(2*sigma^2);
g = exp(-g);
end