close all
clear all

isCompute = 0
isShow = 1;
for ex = 1
    switch ex
        
        case 1
            dName = 'Sergei_2_';
            numFile = 20;
            sigma = 30;
            numAtomOutlier = 0;
            R_low2 = 0; R_high2 = 0;
            R_low3 = 0; R_high3 = 0;
            spExtention = 0.2;
            fqThre = 0;
            extention = 0.2;
            rad = 1;
            threBD = 0.0;
            energyThre = 2.0;
            patchSize = 512;
            t_sc = 1; s_sc = 0.8;
            isCompDeform = 1; red = [10,4];
            fudgeFactor = 0.5; div = 2;
            typeBD = 2; NB = 15;
            isRectangle = 0; isCheckType = 3; isLake = 1;
    end
    
    
    count = 1;
    pic = figure('position',[100 100 850 600]);
    M = moviein(numFile,pic);
    
    for cntFile = 1:numFile
        numAtom = 25;
        fileName = [dName num2str(cntFile) '.png'];
        saveName = [dName num2str(cntFile)];
        fName = sprintf('./results/ex_%s.mat',saveName);
        img = double(imread(fileName));
        if (length(size(img))==3)
            img = img(:,:,1);
        end
        phi = img;
        if isCompute
            
            % denoising
            
            [PSNR, phi] = BM3D(1, phi, sigma);
            
            % Preprocess the image so that we can make the parameters in this
            % code relatively data independent
            phi = phi/max(max(abs(phi)));
            phi = phi - mean(phi(:));
            
            %set up parameters for ss transform
            %N = 1024;  patch size
            %typeBD = 2; NB = [1 num_angle]; where num_angle is the number of grid points in the angle axis
            paraSST = struct('div',div,'NB',[1 NB],'N',patchSize,'red',red,'epsl',1e-4,'is_real',0,'t_sc',t_sc,'s_sc',s_sc,'num_direction',1,...
                'R_low2',R_low2,'R_low3',R_low3,'R_high2',R_high2,'R_high3',R_high3,'SPg',[],'spExtention',spExtention,'extention',extention,...
                'rad',rad,'threBD',threBD,'energyThre',energyThre,'isLake',isLake,'isCheckType',isCheckType,'isRectangle',isRectangle,...
                'isCompDeform',isCompDeform','typeBD',typeBD,'numAtomOutlier',numAtomOutlier);
            
            [orgm,orgn] = size(phi);
            m = orgm; n = orgn;
            
            
            % Estimate basis information, generate filtered images
            [~,info,paraSST] = crystalInfo(phi,paraSST,1,numAtom,extention,fudgeFactor,fqThre*2);
            % info.numPix is the area of numAtom atoms
            %close all;
            results = mainAnalysisSST(phi,paraSST,info,0);
            save(fName,'results','paraSST','info','-v7.3');
        else
            load(fName);
        end
        
        
        
        
        
        spic = subplot(2,2,1);
        imagesc(cutSize2(phi,1/20));axis image;axis off;
        colormap(spic,'gray');
        title('Crystal image');
        set(gca, 'FontSize', 16);
        bb=get(gca);
        set(bb.XLabel, 'FontSize', 16);set(bb.YLabel, 'FontSize', 16);set(bb.ZLabel, 'FontSize', 16);set(bb.Title, 'FontSize', 16);
        
        
        spic = subplot(2,2,2);
        imagesc(cutSize2(results.BD,1/20));
        colormap(spic,1-gray);axis image;axis off;caxis([0 1]);
        title('Defects');
        set(gca, 'FontSize', 16);
        bb=get(gca);
        set(bb.XLabel, 'FontSize', 16);set(bb.YLabel, 'FontSize', 16);set(bb.ZLabel, 'FontSize', 16);set(bb.Title, 'FontSize', 16);
        
        
        pos = find(results.space>0);
        oriMask = results.ori;
        %        oriMask(pos) = 0;
        
        spic = subplot(2,2,3);
        imagesc(cutSize2(oriMask,1/20));axis off;axis image;colorbar;
        colormap(spic,pmkmp());
        title('orientation');
        set(gca, 'FontSize', 16);
        bb=get(gca);
        set(bb.XLabel, 'FontSize', 16);set(bb.YLabel, 'FontSize', 16);set(bb.ZLabel, 'FontSize', 16);set(bb.Title, 'FontSize', 16);
        
        
        spic = subplot(2,2,4);
        imagesc(cutSize2(results.space,1/20));axis off;axis image;caxis([0 1]);
        colormap(spic,'gray');
        title('Vacancy');
        head = sprintf('./results/%s_space.fig',saveName);
        saveas(pic,head);
        head = sprintf('./results/%s_space.png',saveName);
        saveas(pic,head);
        
        
        M(count) = getframe(pic);
        count = count + 1;
    end
    
    
    [h, w, p] = size(M(1).cdata);  % use 1st frame to get dimensions
    hf = figure;
    % resize figure based on frame's w x h, and place at (150, 150)
    set(hf, 'position', [150 150 w h]);
    axis off
    movie(hf,M,1,1);
    movie2avi(M, 'results/Sergei_2.avi', 'compression', 'None','fps',2);
    implay('Sergei_2.avi',2);
end
