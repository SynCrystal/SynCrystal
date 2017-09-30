function [kb avgdx avgdy kk1 kk2] = SS_cswft2_polar_v1(num_direction,fff,h,NB,sigma,R_low,R_high,epsl)
%Input:
%   fff is the image
%   h          h(1) is the spacing of samples in the vertical direction in the image
%              h(2) is the spacing of samples in the horizontal direction in the image
%   NB(1) is the size of radial grid in [R_low R_high]
%   NB(2) is the size of angular grid in [0,pi)
%   sigma       width parameter of the Gaussian window
%   R_low: lower bound of interested frequency
%   R_high: upper bound of interested frequency
%   epsl: the threshold for small coefficients of general curvelet transform
%   Note: The accuracy of the synchrosqueezed transform is determined by the
%   parameter epsl.
%
%Output:
%   kb the synchrosqueezed energy distribution in a polar coordinate
%   size(kb) = NB,
%   cc a tensor storing the windowed Fourier coefficients
%
%   kk1: local wave vector in the first direction
%   kk2: local wave vector in the second direction
%
%by Haizhao Yang

%set up size of image
[Nx Ny] = size(fff);
N = max(Nx,Ny);

%cc is the windowed Fourier coefficients
%aa and bb are the windowed Fourier coefficients with the derivatives in b_1 and b_2
cc = cswft2_fwd(fff,h,sigma);
aa = cswft2_ext_1(fff,h,sigma);
bb = cswft2_ext_2(fff,h,sigma);

EXT = 10^10;
szc = size(cc);
cc = reshape(cc,[szc(1)*szc(2),szc(3),szc(4)]);
aa = reshape(aa,[szc(1)*szc(2),szc(3),szc(4)]);
bb = reshape(bb,[szc(1)*szc(2),szc(3),szc(4)]);
szc = size(cc);
% kk1 local wave vector in the first direction
% kk2 local wave vector in the second direction
gud = find(abs(cc)>epsl);
aa = (real( aa(gud)./cc(gud) / (2*pi*i)));
bb = (real( bb(gud)./cc(gud) / (2*pi*i)));
dist = sqrt(aa.^2+bb.^2);
good = find(dist<=min(N/sqrt(2),R_high) & dist>=R_low);

kk1 = repmat(EXT,szc);
kk1(gud(good)) = aa(good);
kk2 = kk1;
kk2(gud(good)) = bb(good);
clear aa bb;


if(1)
    dr = (R_high-R_low)/(NB(1));
    da = pi/(NB(2))/num_direction;
    %WB = NB/R_high;%2*min(N,2*R_high)/NB;
    kb = zeros(NB(1)+1,NB(2)+1,szc(2),szc(3));
    kb = SS_polar_v1(cc,kk1,kk2,EXT,pi/num_direction,da,dr,NB(1)+1,NB(2)+1,R_low);
end


