function [ss_energy coefTensor LocWavVecx LocWavVecy] = ss_cswft2_fwd(img,h,R_low,R_high,NB,sigma,epsl)
% ss_cswft2_fwd.m - 2D Synchrosqueezed windowed Fourier Transform
%
% Inputs:
%   img         given image of size Nx by Ny
%
% Optional inputs:
%   h          h(1) is the spacing of samples in the vertical direction in the image
%              h(2) is the spacing of samples in the horizontal direction in the image
%   R_low       lower bound of interested frequency, not necessary in this code
%   R_high      upper bound of interested frequency
%   NB          (2*NB+1) by (2*NB+1) is the size of samples in the phase in [-R_high,R_high]*[-R_high,R_high]
%   sigma       width parameter of the Gaussian window
%   epsl        threshold for synchrosqueezed curvelet transform
%               Note: The accuracy of the synchrosqueezed transform is determined by the
%               parameter epsl. However, due to numerical issues, epsl
%               cannot be too small.
%
% Outputs:
%   ss_energy          synchrosqueezed energy distribution
%   coefTensor         a tensor storing the windowed Fourier transform coefficients
%   LocWavVecx         local wave vector in the first direction
%   LocWavVecy         local wave vector in the second direction
%
%by Haizhao Yang


[Nx Ny] = size(img);
N = max(Nx,Ny);
if nargin < 2, h = [1 1]; end;
%if nargin < 3, R_low = 0; end;
R_low = 0;
if nargin < 4, R_high = ceil(N/2); end;
if nargin < 5, NB = ceil(R_high/8); end;
if nargin < 6, sigma = 10; end;
if nargin < 7, epsl = 1e-2; end;

%coefTensor is the windowed Fourier coefficients
%aa and bb are the windowed Fourier coefficients with the derivatives in b_1 and b_2
coefTensor = cswft2_fwd(img,h,sigma);
aa = cswft2_ext_1(img,h,sigma);
bb = cswft2_ext_2(img,h,sigma);

EXT = 10^10;
szc = size(coefTensor);
coefTensor = reshape(coefTensor,[szc(1)*szc(2),szc(3),szc(4)]);
aa = reshape(aa,[szc(1)*szc(2),szc(3),szc(4)]);
bb = reshape(bb,[szc(1)*szc(2),szc(3),szc(4)]);
szc = size(coefTensor);
% LocWavVecx local wave vector in the first direction
% LocWavVecy local wave vector in the second direction
gud = find(abs(coefTensor)>epsl);
aa = (real( aa(gud)./coefTensor(gud) / (2*pi*i)));
bb = (real( bb(gud)./coefTensor(gud) / (2*pi*i)));
dist = sqrt(aa.^2+bb.^2);
good = find(dist<=min(N/sqrt(2),R_high) & dist>=R_low);

LocWavVecx = repmat(EXT,szc);
LocWavVecx(gud(good)) = aa(good);
LocWavVecy = LocWavVecx;
LocWavVecy(gud(good)) = bb(good);
clear aa bb;

if(1)
    WB = NB/R_high;
    ss_energy = zeros(2*NB+1,2*NB+1,szc(2),szc(3));
    for a=1:szc(2)
        for b=1:szc(3)
            tc = coefTensor(a,b,:);
            tk1 = LocWavVecx(a,b,:);
            tk2 = LocWavVecy(a,b,:);
            tc = tc(:);
            gud = find(tk1<EXT);
            tk1 = tk1(gud);
            tk2 = tk2(gud);
            for g=1:length(gud)
                loc1 = round(tk1(g)*WB)+NB+1;
                loc2 = round(tk2(g)*WB)+NB+1;
                ss_energy(loc1,loc2,a,b) = ss_energy(loc1,loc2,a,b) + abs(tc(gud(g))).^2;
            end
        end
    end
end


