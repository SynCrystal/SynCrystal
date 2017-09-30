close all
clear all

%% Set up parameters for the synchrosqueezed curvelet transform
%set up size of image
N = 1024;

%set up an image
%generate the spatial grid
xo = [0:N-1]/N;
fff = zeros(1,N);

amp = 0.0;
F1 = 30;
F2 = 60;
%The local wave number is about 85

xx = xo + amp*sin(2*pi*xo);
f1 = exp(2*pi*i*(F1*xx));

yy = xo + amp*cos(2*pi*xo);
f2 = exp(2*pi*i*(F2*yy));

NM = 0;
%ADD BOTH REAL AND IMAGINARY NOISE
ns = NM*(randn(1,N)+i*randn(1,N));
fff = (f1+f2+ns);

fff = f1;

h = 8;
R_low=0; R_high = N/2; 
NB = R_high/4;
sigma=10;
cc = cswft1_fwd(fff,h,sigma);
aa = cswft1_ext(fff,h,sigma);
figure;plot(abs(cc(:,100)));
figure;plot(abs(aa(:,100)));

epsl=1e-8;
szc = size(cc);
gud = find(abs(cc)>epsl);
aaa = (real( aa(gud)./cc(gud) / (2*pi*i)));

dist = abs(aaa);
good = find(dist<=min(N/sqrt(2),R_high) & dist>=R_low);

instFq = repmat(0,szc);
instFq(gud(good)) = aaa(good);

figure;plot(instFq(:,10));

asdfs
display('begin SSWFT');
if 1
    [ss_energy,~,~,~] = ss_cswft2_fwd(fff,h,R_low,R_high,NB,sigma);
else
    [ss_energy,~,~,~,~] = ss_ct2_fwd(fff);
end

figure;imagesc(abs(cc(:,:,10,10)));
figure;imagesc(ss_energy(:,:,10,10));
