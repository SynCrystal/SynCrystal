function [skeleton,kb,avgdx,avgdy,kk1,kk2,cc,timeSST,timeSkt] = skeletonPolar(num_direction,fff,SPg,NB,rad,is_real,R_low,R_high,epsl,red,t_sc, s_sc,deformFactor,numAgl,isScale,is_cos)
%Input:
%fff is the image
%SPg(1) is the number of samples in the vertical direction in the image
%SPg(2) is the number of samples in the horizontal direction in the image
%NB(1) is the size of radial grid in [R_low R_high]
%NB(2) is the size of angular grid in [0,pi)
%rad is the smallest size of supports of wave packets
%is_real = 0: fff is complex
%is_real = 1: fff is real
%R_low: lower bound of interested frequency
%R_high: upper bound of interested frequency
%epsl: the threshold for small coefficients of general curvelet transform
%Note: The accuracy of the synchrosqueezed transform is determined by the
%parameter epsl.
%   is_cos      Type of the window function
%                   0: C^infinity window function
%                   1: cosine window function
%               [default set to 0]
%   t_sc    scaling parameter for radius
%               [default set to 1-1/8]
%   s_sc    scaling parameter for angle
%               [default set to 1/2+1/8]
%   deformFactor   in (0,1], how much deformation allowed for
%   classification
%
%Output:
%kb the synchrosqueezed energy distribution in a polar coordinate
%   size(kb) = NB,
%ccc the general curvelet transform coefficients
%cc a tensor storing the general curvelet transform coefficients
%
%kk1: local wave vector in the first direction
%kk2: local wave vector in the second direction
%
%by Haizhao Yang

if nargin < 10, red = 1; end;
if nargin < 11, t_sc = 1 - 1/8; end;
if nargin < 12, s_sc = 1/2 + 1/8; end;
if nargin < 13, deformFactor = 0.5; end; 
if nargin < 14, numAgl = 9; end;
if nargin < 15, isScale = 1; end;
if nargin < 16, is_cos = 1; end;

%set up size of image
[Nx Ny] = size(fff);
N = max(Nx,Ny);
%set up size of samples in the space
SPx = SPg(1); SPy = SPg(2);

%ccc is the general curvelet coefficients
%aaa and bbb are the general curvelet coefficients with the derivatives in b_1 and b_2
tic;
[ccc aaa bbb] = gdct2_fwd_red(fff,is_real,[SPx SPy],R_high,R_low,rad,is_cos, t_sc, s_sc, red);

ncl = numel(ccc);
[t1 t2] = size(ccc{1}{1,1});
nclp = 0;
for g = 1:ncl
    nclp = nclp + numel(ccc{g});
end
cc = zeros(t1,t2,nclp);
aa = zeros(t1,t2,nclp);
bb = zeros(t1,t2,nclp);
cnt_nclp = 1;
for g=1:ncl
    [szmccc,sznccc] = size(ccc{g});
    for cnt1 = 1:szmccc
        for cnt2 = 1:sznccc
            cc(:,:,cnt_nclp) = ccc{g}{cnt1,cnt2};
            aa(:,:,cnt_nclp) = aaa{g}{cnt1,cnt2};
            bb(:,:,cnt_nclp) = bbb{g}{cnt1,cnt2};
            cnt_nclp = cnt_nclp + 1;
        end
    end
end
clear aaa bbb ccc;
EXT = 10^10;
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
    % use fine grid to get accurate SST in the polar coordinate
    dr = (R_high-R_low)/(NB(1));
    da = pi/(NB(2))/num_direction;
    %WB = NB/R_high;%2*min(N,2*R_high)/NB;
%     kb = zeros(NB(1)+1,NB(2)+1,size(cc,1),size(cc,2));
%     avgdx = zeros(NB(1)+1,NB(2)+1,size(cc,1),size(cc,2));
%     avgdy = zeros(NB(1)+1,NB(2)+1,size(cc,1),size(cc,2));
%     SS_polar_v2(cc,kk1,kk2,EXT,pi/num_direction,da,dr,NB(1)+1,NB(2)+1,R_low,kb,avgdx,avgdy);
    [kb,avgdx,avgdy] = SS_polar(cc,kk1,kk2,EXT,pi/num_direction,da,dr,NB(1)+1,NB(2)+1,R_low);
    pos = find(kb>0);
    avgdx(pos) = avgdx(pos)./kb(pos);
    avgdy(pos) = avgdy(pos)./kb(pos);
end
timeSST = toc;

if (1) % compute the time-frequency skeleton using a coarse grid
    % find scaleN and shift in angle
    tic;
    [scaleN,shiftAgl] = skeletonPeakShift(kb,avgdx,avgdy);
    pos = find(scaleN~=0);
    minScale = min(scaleN(pos));
    maxScale = max(scaleN(pos));
    if isScale
    % scale for scale invariance
    for cnt1 = 1:szc(1)
        for cnt2 = 1:szc(2)
            if scaleN(cnt1,cnt2)>0
                kk1(cnt1,cnt2,:) = kk1(cnt1,cnt2,:)*minScale/scaleN(cnt1,cnt2);
                kk2(cnt1,cnt2,:) = kk2(cnt1,cnt2,:)*minScale/scaleN(cnt1,cnt2);
            else
                kk1(cnt1,cnt2,:) = EXT;
                kk2(cnt1,cnt2,:) = EXT;
            end
        end
    end
    end
    dr = deformFactor*minScale;
    numRad = ceil(R_high/dr);
    da = pi/numAgl;% 20 degrees
    %skeleton = zeros(numRad,numAgl,t1,t2);
    %     if 1
    skeleton = skeletonPolarMex(cc,kk1,kk2,EXT/maxScale,pi,da,dr,numRad,numAgl,0);
    %     else
    %         skeleton = zeros(numRad,numAgl,szc(1),szc(2));
    %         for ai=1:szc(1)
    %             for bi=1:szc(2)
    %                 for ci=1:szc(3)
    %                     if (kk1(ai,bi,ci)<EXT/maxScale)
    %                         r = sqrt(kk1(ai,bi,ci).^2+kk2(ai,bi,ci).^2);
    %                         if (kk1(ai,bi,ci)>=0)
    %                             agl = mod(acos(kk2(ai,bi,ci)/r),pi);
    %                         else
    %                             agl = mod(3.1415926-acos(kk2(ai,bi,ci)/r),pi);
    %                         end
    %                         loc1 = round(r);
    %                         if (loc1<=0)
    %                             loc1 = loc1 + 1;
    %                         end
    %                         loc2 = ceil(agl/da);
    %                         if (loc2<=0)
    %                             loc2 = loc2 + 1;
    %                         end
    %                         temp_energy = abs(cc(ai,bi,ci));
    %                         skeleton(loc1,loc2,ai,bi) = skeleton(loc1,loc2,ai,bi) + temp_energy;
    %                     end
    %                 end
    %             end
    %         end
    %     end
    
    % shift for rotation invariance
    for cnt1 = 1:szc(1)
        for cnt2 = 1:szc(2)
            shift = max(1,ceil(numAgl*shiftAgl(cnt1,cnt2)/pi))-1;
            skeleton(:,:,cnt1,cnt2) = circshift(skeleton(:,:,cnt1,cnt2),-shift,2);
        end
    end
    timeSkt = toc;
end


