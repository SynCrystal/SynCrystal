close all
clear all

isCompute = 1
isShow = 1;
for ex = [11]
    switch ex
        case 1
            dName = 'ck1.gif';
            R_low2 = 0; R_high2 = 0;
            R_low3 = 0; R_high3 = 0;
            spExtention = 0;
            fqThre = 0;
            extention = 0.0;
            rad = 1;
            threBD = 0.1;
            energyThre = 2.8;
            patchSize = 256;
            t_sc = 1; s_sc = 0.8;
            isCompDeform = 1; red = [4,4];
            fudgeFactor = 0.5; div = 1;
            typeBD = 2; NB = 45;
            isRectangle = 0; isCheckType = 3; isLake = 1;
        case 2
            dName = 'LS1.jpg';
            R_low2 = 0; R_high2 = 0;
            R_low3 = 0; R_high3 = 0;
            spExtention = 0;
            fqThre = 3;
            extention = 0.2;
            rad = 1.6;
            threBD = 0.5;
            energyThre = 2.0;
            patchSize = 512;
            t_sc = 1; s_sc = 1;
            isCompDeform = 1; red = [10,4];
            fudgeFactor = 0.5; div = 2;
            typeBD = 2; NB = 45;
            isRectangle = 0; isCheckType = 3; isLake =1;
        case 3
            dName = 'HVS-Direct-5A.tif';
            R_low2 = 0; R_high2 = 0;
            R_low3 = 0; R_high3 = 0;
            spExtention = 0;
            fqThre = 20;
            extention = 0.4;
            rad = 1.5;
            threBD = 0.5;
            energyThre = 2.5;
            patchSize = 512*2;
            t_sc = 1; s_sc = 0.8;
            isCompDeform = 0; red = [1,1];
            fudgeFactor = 0.5; div = 4;
            typeBD = 2; NB = 45;
            isRectangle = 0; isCheckType = 0; isLake =1;
        case 4
            dName = 'pts3.jpg';
            R_low2 = 0; R_high2 = 0;
            R_low3 = 0; R_high3 = 0;
            spExtention = 0;
            fqThre = 0;
            extention = 0.2;
            rad = 1.5;
            threBD = 0.3;
            energyThre = 2.0;
            patchSize = 512;
            t_sc = 0.8; s_sc = 0.7;
            isCompDeform = 0; red = [4,4];%[8,4] too expensive
            fudgeFactor = 0.5; div = 4;
            typeBD = 2; NB = 45;
            isRectangle = 0; isCheckType = 0; isLake =0;
        case 5
            dName = 'ITF2.png';
            R_low2 = 0; R_high2 = 0;
            R_low3 = 0; R_high3 = 0;
            spExtention = 0;
            fqThre = 0;
            extention = 0.5;
            rad = 1.5;
            threBD = 0.0;
            energyThre = 2.0;
            patchSize = 512;
            t_sc = 1; s_sc = 0.8;
            isCompDeform = 0; red = [8,4];
            fudgeFactor = 0.5; div = 2;
            typeBD = 1; NB = 45;
            isRectangle = 0; isCheckType = 0; isLake =0;
            
        case 6
            dName = '0000731.tif';
            R_low2 = 0; R_high2 = 0;
            R_low3 = 0; R_high3 = 0;
            spExtention = 0;
            fqThre = 0;
            extention = 0.0;
            rad = 0.8;
            threBD = 0.0;
            energyThre = 2.0;
            patchSize = 512;
            t_sc = 0.8; s_sc = 0.65;
            isCompDeform = 1; red = [5,1];
            fudgeFactor = 0.5; div = 2;
            typeBD = 2; NB = 45;
            isRectangle = 0; isCheckType = 0; isLake =0;
        case 7
            dName = '0001025.tif';
            R_low2 = 0; R_high2 = 0;
            R_low3 = 0; R_high3 = 0;
            spExtention = 0;
            fqThre = 0;
            extention = 0.0;
            rad = 1;
            threBD = 0.0;
            energyThre = 2.0;
            patchSize = 512;
            t_sc = 1; s_sc = 0.8;
            isCompDeform = 1; red = [10,1];
            fudgeFactor = 0.5; div = 2;
            typeBD = 2; NB = 45;
            isRectangle = 0; isCheckType = 0; isLake =0;
        case 8
            dName = '0005281.tif';
            R_low2 = 0; R_high2 = 0;
            R_low3 = 0; R_high3 = 0;
            spExtention = 0;
            fqThre = 0;
            extention = 0.0;
            rad = 1;
            threBD = 0.0;
            energyThre = 2.0;
            patchSize = 512;
            t_sc = 1; s_sc = 0.8;
            isCompDeform = 1; red = [10,1];
            fudgeFactor = 0.5; div = 2;
            typeBD = 2; NB = 45;
            isRectangle = 0; isCheckType = 0; isLake =0;
        case 9
            dName = '0013786.tif';
            R_low2 = 0; R_high2 = 0;
            R_low3 = 0; R_high3 = 0;
            spExtention = 0;
            fqThre = 0;
            extention = 0.0;
            rad = 1;
            threBD = 0.0;
            energyThre = 2.0;
            patchSize = 512;
            t_sc = 1; s_sc = 0.8;
            isCompDeform = 1; red = [10,1];
            fudgeFactor = 0.5; div = 2;
            typeBD = 2; NB = 45;
            isRectangle = 0; isCheckType = 0; isLake =0;
        case 10
            dName = '0013842.tif';
            R_low2 = 0; R_high2 = 0;
            R_low3 = 0; R_high3 = 0;
            spExtention = 0;
            fqThre = 0;
            extention = 0.0;
            rad = 1;
            threBD = 0.0;
            energyThre = 2.0;
            patchSize = 512;
            t_sc = 1; s_sc = 0.8;
            isCompDeform = 1; red = [10,1];
            fudgeFactor = 0.5; div = 2;
            typeBD = 2; NB = 45;
            isRectangle = 0; isCheckType = 0; isLake =0;
        case 11
            dName = '0013891.tif';
            R_low2 = 0; R_high2 = 0;
            R_low3 = 0; R_high3 = 0;
            spExtention = 0;
            fqThre = 0;
            extention = 0.1;
            rad = 1.6;
            threBD = 0.0;
            energyThre = 2.0;
            patchSize = 512*2;
            t_sc = 1; s_sc = 0.8;
            isCompDeform = 1; red = [8,4];%[10,4];
            fudgeFactor = 0.5; div = 2;
            typeBD = 2; NB = 45;
            isRectangle = 0; isCheckType = 3; isLake =0
    end
    
    img = double(imread(dName));
    if (length(size(img))==3)
        img = img(:,:,1);
    end
    dName
    
    % Preprocess the image so that we can make the parameters in this
    % code relatively data independent
    phi = img;%(:,end-235:end);
    phi = phi/max(max(abs(phi)));
    phi = phi - mean(phi(:));
    
    %set up parameters for ss transform
    %N = 1024;  patch size
    %typeBD = 2; NB = [1 num_angle]; where num_angle is the number of grid points in the angle axis
    paraSST = struct('div',div,'NB',[1 NB],'N',patchSize,'red',red,'epsl',1e-4,'is_real',0,'t_sc',t_sc,'s_sc',s_sc,'num_direction',1,...
        'R_low2',R_low2,'R_low3',R_low3,'R_high2',R_high2,'R_high3',R_high3,'SPg',[],'spExtention',spExtention,'extention',extention,...
        'rad',rad,'threBD',threBD,'energyThre',energyThre,'isLake',isLake,'isCheckType',isCheckType,'isRectangle',isRectangle,...
        'isCompDeform',isCompDeform','typeBD',typeBD);
    
    [orgm,orgn] = size(phi);
    m = orgm; n = orgn;
    
    numAtom = 25;
    fName = sprintf('./results/ex%d.mat',ex);
    if isCompute
        % Estimate basis information, generate filtered images
        [~,info,paraSST] = crystalInfo(phi,paraSST,1,numAtom,extention,fudgeFactor,fqThre*2);
        % info.numPix is the area of numAtom atoms
        %close all;
        results = mainAnalysisSST(phi,paraSST,info,0);
        save(fName,'results','paraSST','info','-v7.3');
    else
        load(fName);
    end
    
    % cut
    if 0
        if isCutm == 1
            results.BD = cutSize(results.BD,1,orgm/paraSST.N);
            results.ori = cutSize(results.ori,1,orgm/paraSST.N);
            results.mask = cutSize(results.mask,1,orgm/paraSST.N);
            results.matType = cutSize(results.matType,1,orgm/paraSST.N);
            results.matLake = cutSize(results.matLake,1,orgm/paraSST.N);
            results.outlier = cutSize(results.outlier,1,orgm/paraSST.N);
            results.segmentation = cutSize(results.segmentation,1,orgm/paraSST.N);
        end
        if isCutn == 1
            results.BD = cutSize(results.BD,2,orgn/paraSST.N);
            results.ori = cutSize(results.ori,2,orgn/paraSST.N);
            results.mask = cutSize(results.mask,2,orgn/paraSST.N);
            results.matType = cutSize(results.matType,2,orgn/paraSST.N);
            results.matLake = cutSize(results.matLake,2,orgn/paraSST.N);
            results.outlier = cutSize(results.outlier,2,orgn/paraSST.N);
            results.segmentation = cutSize(results.segmentation,2,orgn/paraSST.N);
        end
    end
    
    if isShow
        pic = figure('Name','Crystal image');imagesc(cutSize2(phi,1/20));axis image;axis off;colormap gray;
        head = sprintf('./results/%d_img.fig',ex);
        saveas(pic,head);
        head = sprintf('./results/%d_img.png',ex);
        saveas(pic,head);
        
        pic = figure('Name','types');imagesc(cutSize2(results.matType,1/20));axis image;axis off;caxis([1 3]);colorbar;
        head = sprintf('./results/%d_type.fig',ex);
        saveas(pic,head);
        head = sprintf('./results/%d_type.png',ex);
        saveas(pic,head);
        
        pic = figure('Name','Grain boundary');imagesc(cutSize2(results.BD,1/20));colormap (1-gray);axis image;axis off;caxis([0 1]);
        head = sprintf('./results/%d_BD.fig',ex);
        saveas(pic,head);
        head = sprintf('./results/%d_BD.png',ex);
        saveas(pic,head);
        
        pos = find(results.matLake>0);
        oriMask = results.ori;
        if 0%isLake
            oriMask(pos) = 0;
            posType = find(results.matType==2);
            oriMask(posType) = mod(oriMask(posType),90);
            posType = find(results.matType==3);
            oriMask(posType) = mod(oriMask(posType),60);
        end
        
        pic = figure('Name','orientation');
        imagesc(cutSize2(oriMask,1/20));axis off;axis image;colorbar;colormap(pmkmp());
        % if isCheckType == 2
        %     caxis([0,90]);
        % else
        %     caxis([0,60]);
        % end
        head = sprintf('./results/%d_ori_mask.fig',ex);
        saveas(pic,head);
        head = sprintf('./results/%d_ori_mask.png',ex);
        saveas(pic,head);
        
        if ex == 17
            pic = figure('Name','orientation');
            imagesc(cutSize2(oriMask(end/4+1:3*end/4,end/4+1:3*end/4),1/20));axis off;axis image;colorbar;colormap(pmkmp());
            head = sprintf('./results/%d_ori_maskz.fig',ex);
            saveas(pic,head);
            head = sprintf('./results/%d_ori_maskz.png',ex);
            saveas(pic,head);
            
            pic = figure('Name','Grain boundary');imagesc(cutSize2(results.BD(end/4+1:3*end/4,end/4+1:3*end/4),1/20));colormap (1-gray);axis image;axis off;caxis([0 1]);
            head = sprintf('./results/%d_BDz.fig',ex);
            saveas(pic,head);
            head = sprintf('./results/%d_BDz.png',ex);
            saveas(pic,head);
            
            pic = figure('Name','Crystal image');imagesc(cutSize2(phi(end/4+1:3*end/4,end/4+1:3*end/4),1/20));axis image;axis off;colormap gray;
            head = sprintf('./results/%d_imgz.fig',ex);
            saveas(pic,head);
            head = sprintf('./results/%d_imgz.png',ex);
            saveas(pic,head);
        end
        
        pic = figure('Name','lakes');imagesc(cutSize2(results.matLake,1/20));axis off;axis image;caxis([0 1]);colorbar;
        head = sprintf('./results/%d_matLake.fig',ex);
        saveas(pic,head);
        head = sprintf('./results/%d_matLake.png',ex);
        saveas(pic,head);
        
        pic = figure('Name','outliers');imagesc(cutSize2(results.outlier,1/20));axis off;axis image;caxis([0 1]);colorbar;
        head = sprintf('./results/%d_outlier.fig',ex);
        saveas(pic,head);
        head = sprintf('./results/%d_outlier.png',ex);
        saveas(pic,head);
        
        pic = figure('Name','grain segmentation');imagesc(cutSize2(results.segmentation,1/20));axis off;axis image;colorbar;colormap(pmkmp());
        head = sprintf('./results/%d_segmentation.fig',ex);
        saveas(pic,head);
        head = sprintf('./results/%d_segmentation.png',ex);
        saveas(pic,head);
        
        if paraSST.isCompDeform
            % volume distortion
            vol = results.deformG(:,:,1).*results.deformG(:,:,4)-results.deformG(:,:,2).*results.deformG(:,:,3);
            %vol(pos) = 1;
            pic = figure('Name','Distortion Volume');imagesc(cutSize2(vol,1/20)-1); axis off; axis image; colorbar;colormap(pmkmp());
            head = sprintf('./results/%d_vol.fig',ex);
            saveas(pic,head);
            head = sprintf('./results/%d_vol.png',ex);
            saveas(pic,head);
            
            [mG nG] = size(results.deformG(:,:,1,1));
            %diff = zeros(mG,nG);
            for cnt = 1:mG
                for cnt2 = 1:nG
                    temp = [results.deformG(cnt,cnt2,1,1) results.deformG(cnt,cnt2,1,2);results.deformG(cnt,cnt2,2,1) results.deformG(cnt,cnt2,2,2)];
                    [U H] = polarDec(temp);
                    [ev ew] = eig(H);
                    diff(cnt,cnt2) = abs(ew(1,1)-ew(2,2));
                end
            end
            
            %diff(posBD) = 0;
            pic = figure('Name','Difference in principle stretches');
            %diff(pos) = 0;
            imagesc(diff); axis off; axis image;set(gca,'xtick',[]);set(gca,'ytick',[]);colorbar;
            %alphamask(bwMask,[0.3 0.3 0.3],1);
            head = sprintf('./results/%d_dist.fig',ex);
            saveas(pic,head);
            head = sprintf('./results/%d_dist.png',ex);
            saveas(pic,head);
        end
    end
end
