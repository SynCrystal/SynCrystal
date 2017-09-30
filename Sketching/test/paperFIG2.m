example = 1; IsFC = 0;


switch example
    case 1 % use high subsampleRate or small red
        load 'GB1.mat';
        if 1
            % less memory than before
            parameters = struct('extendRatio',0.0,'threRatio',0.1, ...
                'epsFreqBand',1e-2,'is_real',0,'numWave',3, ...
                'numGridPerSector',30,'red',[4,1],'rad',2, ...
                's_sc',0.95,'t_sc',0.95,'epsSST',1e-1);
            subsampleRate = 4;
        else
            % much less memory than before
            parameters = struct('extendRatio',0.0,'threRatio',0.1, ...
                'epsFreqBand',1e-2,'is_real',0,'numWave',3, ...
                'numGridPerSector',30,'red',[4,1],'rad',2, ...
                's_sc',0.75,'t_sc',0.95,'epsSST',1e-1);
            subsampleRate = 4;
        end
        energyThre = 2.3;
        disVolThre = 0.1;
        numIter = 40;
        if ~IsFC
            coarseThre = 0.3; fineThre = 0.3;  %we can also choose coarseThre = fineThre = 0.5 or 0.3, larger defect region
            fprintf('Using coarseThre = fineThre = %f\n',fineThre);
        else
            coarseThre = 0.2; fineThre = 0.3;
            fprintf('Using coarseThre %f;   fineThre = %f\n',coarseThre,fineThre);
        end
        
end

% normalize the given image so as to use uniform thresholding parameters
fff = phi;
fff = fff/max(abs(fff(:)));
fff = fff - mean(fff(:));

%% -------------------------------------------------------------------------
% The synchrosqueezed wave packet transform provids an initial guess
tic;
plotFIG = false;
[radii angles TTEng_1st TTEng_2nd masses ss_energy] = initialGuess(fff,subsampleRate,plotFIG,parameters); % assumes y-coordinate pointing into direction of growing first index

% bring wave vectors into order according to angle
[angles,perm] = sort(angles,3);
for k = 1:size(perm,1)
    for l = 1:size(perm,2)
        radii(k,l,:) = radii(k,l,perm(k,l,:));
    end
end

%% Show angles and radii
if 0
    figure;
    subplot(2,3,1);imagesc(angles(:,:,1));axis image;colorbar; title('new angle 1');
    subplot(2,3,2);imagesc(angles(:,:,2));axis image;colorbar; title('new angle 2');
    subplot(2,3,3);imagesc(angles(:,:,3));axis image;colorbar; title('new angle 3');
    subplot(2,3,4);imagesc(radii(:,:,1));axis image;colorbar; title('new radii 1');
    subplot(2,3,5);imagesc(radii(:,:,2));axis image;colorbar; title('new radii 2');
    subplot(2,3,6);imagesc(radii(:,:,3));axis image;colorbar; title('new radii 3');
end

% find defect region and their connected components
if 1
    binaryIndexCoarse = 1-findDeftArea(smoothImage( sum(masses,3), 4, 1),energyThre,coarseThre);
    [CCCoarse,colComponentsCoarse] = periodicConnComp( binaryIndexCoarse );
    binaryIndex = 1-findDeftArea(smoothImage( sum(masses,3), 4, 1),energyThre,fineThre);
end

[CC,colComponents] = periodicConnComp( binaryIndex );

if 1
    massfig = smoothImage( sum(masses,3), 4, 1);
    [szm szn ] = size(massfig);
    temp = massfig(round(391/1024*szm):round(518/1024*szm),round(470/1024*szn):round(598/1024*szn));
    pic = figure;imagesc(temp);colorbar;axis image;axis off;colormap(brewermap([],'*RdYlBu'));
    saveas(pic,'./results/patchmass.fig');
    
        pos1 = [round((round(391/1024*szm)+round(518/1024*szm))/2),round((round(470/1024*szn)+round(598/1024*szn))/2)];
        pic = figure;imagesc([0 180],[148 204],ss_energy(:,:,pos1(1),pos1(2)));axis square;colorbar; colormap (1-gray);
        xlabel('Angle');ylabel('Wave Number');axis xy;set(gca,'xgrid','on')
        set(gca,'xtick',[0:60:180]);
        saveas(pic,'./results/patchSS1.fig');
    
        pos1 = [round((round(391/1024*szm)+round(518/1024*szm))/2),round((round(470/1024*szn)+round(598/1024*szn))/2)-40];
        pic = figure;imagesc([0 180],[148 204],ss_energy(:,:,pos1(1),pos1(2)));axis square;colorbar; colormap (1-gray);
        xlabel('Angle');ylabel('Wave Number');axis xy;set(gca,'xgrid','on')
        set(gca,'xtick',[0:60:180]);
        saveas(pic,'./results/patchSS2.fig');

        temp = binaryIndex(round(391/1024*szm):round(518/1024*szm),round(470/1024*szn):round(598/1024*szn));
        pic = figure;imagesc(temp);colorbar;axis image;axis off;colormap gray;
        saveas(pic,'./results/patchbw.fig');
end
toc
%% -------------------------------------------------------------------------
tic;
% find initial noisy deformation gradient, assuming y-coordinate pointing into direction of decreasing 1st index, and using two different point group interpretations
G1 = initGFromSSPeaks( angles, radii, pi/3*[0;1;2], median(radii(:))*ones(3,1) );

%% -------------------------------------------------------------------------
% find locations of point group effect (first define the point group)
rot = [cos(pi/3) sin(pi/3);-sin(pi/3) cos(pi/3)];
pointGroup = eye(2);
for j = 1:5
    pointGroup = cat(3,pointGroup,pointGroup(:,:,j)*rot);
end
[L,G0] = findPointGroupJumpsBySweeping( pointGroup, G1, colComponentsCoarse );

address = sprintf('./results/');

[m n] = size(G0(:,:,1,1));
bwMask = 1 - binaryIndex;
posBD = find(bwMask==1);
agl = zeros(m,n);
diff = zeros(m,n);
for cnt = 1:m
    for cnt2 = 1:n
        temp = [G0(cnt,cnt2,1,1) G0(cnt,cnt2,1,2);G0(cnt,cnt2,2,1) G0(cnt,cnt2,2,2)];
        [U H] = polarDec(temp);
        agl(cnt,cnt2) = mod(acos(U(1,1)),pi/3)*180/pi;
        [ev ew] = eig(H);
        diff(cnt,cnt2) = abs(ew(1,1)-ew(2,2));
    end
end
agl(posBD) = 0;
pic = figure;imagesc(agl(round(391/1024*szm):round(518/1024*szm),round(470/1024*szn):round(598/1024*szn))); axis image;set(gca,'xtick',[]);set(gca,'ytick',[]);
colorbar;colormap(pmkmp());
alphamask(bwMask(round(391/1024*szm):round(518/1024*szm),round(470/1024*szn):round(598/1024*szn)),[0.3 0.3 0.3],1);
name = sprintf('%sAglPatch.fig',address);
saveas(pic,name);

diff(posBD) = 0;
pic = figure;imagesc(diff(round(391/1024*szm):round(518/1024*szm),round(470/1024*szn):round(598/1024*szn))); axis image;set(gca,'xtick',[]);set(gca,'ytick',[]);
colorbar;colormap(brewermap([],'*RdYlBu'));
alphamask(bwMask(round(391/1024*szm):round(518/1024*szm),round(470/1024*szn):round(598/1024*szn)),[0.3 0.3 0.3],1);
name = sprintf('%sDiffPatch.fig',address);
caxis([0 0.1]);
saveas(pic,name);


pic = figure;
% volume distortion
vol = G0(:,:,1).*G0(:,:,4)-G0(:,:,2).*G0(:,:,3);
imagesc(vol(round(391/1024*szm):round(518/1024*szm),round(470/1024*szn):round(598/1024*szn))-1); axis image;  h=colorbar; caxis([-0.1 0.1]);
colorbar;colormap(brewermap([],'*RdYlBu'));
title('Distortion Volume');
set(gca,'xtick',[]);set(gca,'ytick',[]);
name = sprintf('%sVolPatch.fig',address);
alphamask(bwMask(round(391/1024*szm):round(518/1024*szm),round(470/1024*szn):round(598/1024*szn)),[0.3 0.3 0.3],1);
saveas(pic,name);


temp = G0(round(391/1024*szm):round(518/1024*szm),round(470/1024*szn):round(598/1024*szn),1);
pic = figure;imagesc(temp);colorbar;axis image;axis off;caxis([-1 1]);colormap(brewermap([],'*RdYlBu'));
saveas(pic,'./results/patchG1.fig');
temp = G0(round(391/1024*szm):round(518/1024*szm),round(470/1024*szn):round(598/1024*szn),2);
pic = figure;imagesc(temp);colorbar;axis image;axis off;caxis([-1 1]);colormap(brewermap([],'*RdYlBu'));
saveas(pic,'./results/patchG2.fig');
temp = G0(round(391/1024*szm):round(518/1024*szm),round(470/1024*szn):round(598/1024*szn),3);
pic = figure;imagesc(temp);colorbar;axis image;axis off;caxis([-1 1]);colormap(brewermap([],'*RdYlBu'));
saveas(pic,'./results/patchG3.fig');
temp = G0(round(391/1024*szm):round(518/1024*szm),round(470/1024*szn):round(598/1024*szn),4);
pic = figure;imagesc(temp);colorbar;axis image;axis off;caxis([-1 1]);colormap(brewermap([],'*RdYlBu'));
saveas(pic,'./results/patchG4.fig');
