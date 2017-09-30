function [ss_energy coefMat instFq LocWavVecy] = ss_cswft2_fwd(f,h,R_low,R_high,NB,sigma,epsl)
% ss_cswft1_fwd.m - 1D Synchrosqueezed windowed Fourier Transform
%
% Inputs:
%   f         given vector of size N
%
% Optional inputs:
%   h          h is the spacing of samples
%   R_low       lower bound of interested frequency, not necessary in this code
%   R_high      upper bound of interested frequency
%   NB          (2*NB+1) is the size of samples in the frequency in [-R_high,R_high]
%   sigma       width parameter of the Gaussian window
%   epsl        threshold for synchrosqueezed curvelet transform
%               Note: The accuracy of the synchrosqueezed transform is determined by the
%               parameter epsl. However, due to numerical issues, epsl
%               cannot be too small.
%
% Outputs:
%   ss_energy          synchrosqueezed energy distribution
%   coefMat            a matrix storing the windowed Fourier transform coefficients
%   instFq             instantaneous frequency estimates
%
%by Haizhao Yang


N = length(f);
if nargin < 2, h = 1; end;
%if nargin < 3, R_low = 0; end;
R_low = 0;
if nargin < 4, R_high = ceil(N/2); end;
if nargin < 5, NB = ceil(R_high/8); end;
if nargin < 6, sigma = 10; end;
if nargin < 7, epsl = 1e-2; end;

%coefMat is the windowed Fourier coefficients
%aa is the windowed Fourier coefficients with the derivatives in time
coefMat = cswft1_fwd(f,h,sigma);
aa = cswft1_ext(f,h,sigma);

EXT = 10^10;
szc = size(coefMat);
% instFq local wave vector in the first direction
gud = find(abs(coefMat)>epsl);
aa = (real( aa(gud)./coefMat(gud) / (2*pi*i)));
dist = abs(aa);
good = find(dist<=min(N/sqrt(2),R_high) & dist>=R_low);

instFq = repmat(EXT,szc);
instFq(gud(good)) = aa(good);
clear aa;

if(1)
    WB = NB/R_high;
    ss_energy = zeros(2*NB+1,szc(2));
    for a=1:szc(2)
        tc = coefMat(a,:);
        tk1 = instFq(a,:);
        gud = find(tk1<EXT);
        tk1 = tk1(gud);
        for g=1:length(gud)
            loc1 = round(tk1(g)*WB)+NB+1;
            ss_energy(loc1,a) = ss_energy(loc1,a) + abs(tc(gud(g))).^2;
        end
    end
end


