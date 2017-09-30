function [phiHighPass,crystalInfo,paraSST] = crystalInfo(phi,paraSST,isShow,numAtom,extention,fudgeFactor,highPassFq,lowPassFq)
% phi: given crystal image
% highPassFq: frequency parameter for high pass filter
% lowPassFq: frequency parameter for low pass filter
% paraSST.N: patch size
% isShow: is showing results
% Main parameter for detection
% fudgeFactor = .5;
if nargin < 8, lowPassFq = 0; end;

%% Compute high pass image for flat field correction
[szPhi1 szPhi2] = size(phi);
phih = fftshift(fft2(phi));
fqThreVtch = round(highPassFq*szPhi1/paraSST.N);
fqThreHrzh = round(highPassFq*szPhi2/paraSST.N);
temp = phih;
temp(round(end/2)-fqThreVtch:round(end/2)+fqThreVtch,round(end/2)-fqThreHrzh:round(end/2)+fqThreHrzh) = 0;
phiHighPass = real(ifft2(ifftshift(temp)));
if isShow
    figure;imagesc(phiHighPass);axis image;title('high pass image');
end

%% Compute low pass image for denoising
if lowPassFq > 0
    fqThreVtcl = round(lowPassFq*szPhi1/paraSST.N);
    fqThreHrzl = round(lowPassFq*szPhi2/paraSST.N);
else
    fqThreVtcl = round(szPhi1/4);
    fqThreHrzl = round(szPhi2/4);
end
temp2 = zeros(szPhi1,szPhi2);
temp2(round(end/2)-fqThreVtcl:round(end/2)+fqThreVtcl,round(end/2)-fqThreHrzl:round(end/2)+fqThreHrzl) = temp(round(end/2)-fqThreVtcl:round(end/2)+fqThreVtcl,round(end/2)-fqThreHrzl:round(end/2)+fqThreHrzl);
phiLowPass = real(ifft2(ifftshift(temp2)));
if isShow
    figure;imagesc(phiLowPass);axis image;title('low pass image');
end

%% Identify space between grain by variation
% Cannot identify the space if there is noise
% Can identify space if there is flat field
% First step with high pass image
dh = phiHighPass-phiHighPass(:,[2:end 1]);
dv = phiHighPass-phiHighPass([2:end 1],:);
dphi = sqrt(dh.^2+dv.^2);
dphi = smoothImage(dphi,3,1);
[~, thresGlobal] = edge(dphi, 'sobel');
pos = find(dphi>thresGlobal*fudgeFactor);
dphi1 = zeros(size(dphi));
dphi1(pos) = 1;
dphi1 = smoothImage(1-dphi1,10,3);
pos = find(dphi1>1/10);
dphi1 = zeros(size(dphi1));
dphi1(pos) = 1;
dphi1 = 1-dphi1;

% Second step with low pass image
dh = phiLowPass-phiLowPass(:,[2:end 1]);
dv = phiLowPass-phiLowPass([2:end 1],:);
dphi = sqrt(dh.^2+dv.^2);
dphi = smoothImage(dphi,3,1);
[~, thresGlobal] = edge(dphi, 'sobel');
pos = find(dphi>thresGlobal*fudgeFactor);
dphi2 = zeros(size(dphi));
dphi2(pos) = 1;
% dphi2 = smoothImage(1-dphi2,10,3);
% pos = find(dphi2>1/10);
% dphi2 = zeros(size(dphi2));
% dphi2(pos) = 1;
% dphi2 = 1-dphi2;
dphi = min(dphi1,dphi2); % use min to combine results

%% Identify space by edge detection with estimated threshold
% Cannot identify space if there is flat field
% Can identify the space if there is noise
% Hence we use high-pass filtered data
[~, threshold] = edge(phiHighPass, 'sobel');
BWs = edge(phiHighPass,'sobel', threshold * fudgeFactor);
se90 = strel('line', 6, 90);
se0 = strel('line', 6, 0);
matLake = imdilate(BWs, [se90 se0]);

% Identify a patch full of grains
space = min(matLake,dphi);
space = smoothImage(space,paraSST.N,paraSST.N/2);
val = max(space(:));
[pos1 pos2] = find(space==val(1));
range1 = (pos1(1)-paraSST.N/2):(pos1(1)+paraSST.N/2-1);
range2 = (pos2(1)-paraSST.N/2):(pos2(1)+paraSST.N/2-1);

phi = phiHighPass(mod(range1-1,szPhi1)+1,mod(range2-1,szPhi2)+1);

% Zero padding for non-periodic image
dm = round(paraSST.N*0.1/paraSST.div)*paraSST.div;
dn = round(paraSST.N*0.1/paraSST.div)*paraSST.div;
phi = [zeros(dm,dn*2+paraSST.N);zeros(paraSST.N,dn) phi zeros(paraSST.N,dn);zeros(dm,dn*2+paraSST.N)];
phi = -phi;
%estimate Frequency band
[m n] = size(phi);
if m ~= n
    mm = max(m,n);
    if m > n
        phi = [phi zeros(m,m-n)];
    else
        phi = [phi; zeros(n-m,n)];
    end
    idx1 = 1:mm;
    idx2 = 1:mm;
    phi = phi(idx1,idx2);
else
    mm = m;
    idx1 = 1:m;
    idx2 = 1:n;
    phi = phi(idx1,idx2);
end
Y=fftshift(fft2(ifftshift(phi)));
if isShow
    pic = figure;imagesc([-mm/2 mm/2],[-mm/2 mm/2],abs(Y));axis square; colormap (1-gray);
end

idx1 = -(m+mod(m+1,2)-1)/2:(m-mod(m+1,2)-1)/2;
idx2 = -(n+mod(n+1,2)-1)/2:(n-mod(n+1,2)-1)/2;
[phi1 phi2] = ndgrid(idx1,idx2);
R = round(sqrt(phi1.^2+phi2.^2));
R = R(:);
Y = Y(:);
thre = max(Y)*1e-4;
pos = find(Y>thre);
L = max(R);
spec = zeros(1,L+1);
for cnt = 1:length(pos)
    spec(R(pos(cnt))+1) = spec(R(pos(cnt))+1) + Y(pos(cnt))*Y(pos(cnt))';
end
spec(1:highPassFq) = 0;

if isShow
    pic = figure;hold on; plot(spec);
end
[val pos] = max(spec);
numPix = m/pos;  %number of pixels for the radius of one atom
numPixAtom = numPix;
numPix = numAtom*numPix^2;% area for numAtom atoms
thre = val*0.1;
ed = pos+1;
cnt = 0;

while cnt == 0 & ed < L
    if spec(ed)>=spec(ed-1) & spec(ed) < thre
        cnt = 1;
    end
    ed = ed + 1;
end
ed = min(L,ceil((ed-1)*(1+extention)));
st = pos -1;
cnt = 0;
while cnt == 0 & st > 1
    if spec(st)>=spec(st+1) & spec(st) < thre
        cnt = 1;
    end
    st = st - 1;
end
st = max(1,floor((st+1)*(1-extention)));
if st <= 2
    spec(1:ed) = 0;
    [val pos] = max(spec);
    thre = val*0.1;
    ed = pos+1;
    cnt = 0;
    while cnt == 0 & ed < L
        if spec(ed)>=spec(ed-1) & spec(ed) < thre
            cnt = 1;
        end
        ed = ed + 1;
    end
    ed = min(L,ceil((ed-1)*(1+extention)));
    st = pos -1;
    cnt = 0;
    while cnt == 0 & st > 1
        if spec(st)>=spec(st+1) & spec(st) < thre
            cnt = 1;
        end
        st = st - 1;
    end
    st = max(1,floor((st+1)*(1-extention)));
end

L = round(2*length(spec)/3);
ed = min(ed,L);
spec = spec(1:L);
if isShow
    plot(st,spec(st),'ro'); plot(ed,spec(ed),'ro');hold off;
    title('Fourier energy distribution and identified energy band');axis square;xlabel('Frequency');ylabel('Energy');
end
%set up size of samples in the space
if (paraSST.R_high2==0 | paraSST.R_high3 ==0)
    paraSST.R_low2 =st; paraSST.R_high2 = ed;
    paraSST.R_low3 =st; paraSST.R_high3 = ed;
end
paraSST.SPg = round([m n]/paraSST.div/4)*4;


% Correct small piece of lakes
[CC,colComponents] = periodicConnComp( dphi1, 0 );
for cnt = 1:CC.NumObjects
    if numel(CC.PixelIdxList{cnt}) < numPixAtom^2
        dphi1(CC.PixelIdxList{cnt}) = 1;
    end
end
% Correct small piece of lakes
[CC,colComponents] = periodicConnComp( dphi2, 0 );
for cnt = 1:CC.NumObjects
    if numel(CC.PixelIdxList{cnt}) < numPixAtom^2
        dphi2(CC.PixelIdxList{cnt}) = 1;
    end
end
dphi = min(dphi1,dphi2); % use min to combine results
% Correct small piece of lakes
[CC,colComponents] = periodicConnComp(matLake, 0 );
for cnt = 1:CC.NumObjects
    if numel(CC.PixelIdxList{cnt}) < numPixAtom^2
        matLake(CC.PixelIdxList{cnt}) = 1;
    end
end
space = min(matLake,dphi);
if isShow
    figure;imagesc(dphi1);colorbar;title('dphi1');
    figure;imagesc(dphi2);colorbar;title('dphi2');
    figure;imagesc(matLake);colorbar;title('matLake');
    figure;imagesc(space);colorbar;title('space');
end

%% Summarize results
crystalInfo = struct('dphi',dphi,'dphi1',dphi1,'dphi2',dphi2,'matLake',matLake,'numPix',numPix,'numPixAtom',numPixAtom,'space',space);
