close all
clear all
warning off;

% This code is to test the classification algorithm by phase space
% sketching using synthetic examples in the paper.

is_show = 0;
R_low2 = 0; R_high2 = 0;
R_low3 = 0; R_high3 = 0;
spExtention = 0;
fqThre = 0;
extention = 0.7;
rad = 1.5;
threBD = 0.3;
energyThre = 2.0;
t_sc = 0.8; s_sc = 0.7;
red = 8;%20
div = 4;
NB = [30,35];
isLake = 0;
gdSzType = 5;
patchSize = 512;


for ex = 1:5
    switch ex
        case 1
            sigma = 0.2; % parameter for spectral clustering
            deformFactor = 0.25;
            distType = 4;
            numAgl = 12;
            clsThre = 0.5;
            isScale = 1;
            
            N = 512;
            x=0:1/N:(N-1)/N;
            [xo yo] = ndgrid(x);
            amp = 0.0;
            xx = xo + amp*sin(2*pi*xo);
            yy = yo + amp*sin(2*pi*yo);
            bmp = 0.0;
            alpha = -pi/8 + bmp*sin(2*pi*xo) + bmp*cos(2*pi*yo);
            F = 50;
            phi = exp(2*pi*i*F*(cos(alpha).*xx+sin(alpha).*yy)) + exp(2*pi*i*F*(cos(alpha+pi/3).*xx+sin(alpha+pi/3).*yy)) + exp(2*pi*i*F*(cos(alpha+2*pi/3).*xx+sin(alpha+2*pi/3).*yy));
            phi = real(phi)/3;
            % pause;
            alpha = -pi/8  + bmp*cos(2*pi*xo) + bmp*sin(2*pi*yo);
            F = 50;
            phi2 = exp(2*pi*i*F*(cos(alpha).*xx+sin(alpha).*yy)) + exp(2*pi*i*F*(cos(alpha+pi/2).*xx+sin(alpha+pi/2).*yy));
            phi2 = (real(phi2)/2);
            phi(end/2+1:end,1:end/2) = phi2(end/2+1:end,1:end/2);
            F = 60;
            phi3 = exp(2*pi*i*F*(cos(alpha).*xx+sin(alpha).*yy)) + exp(2*pi*i*F*(cos(alpha+pi/4).*xx+sin(alpha+pi/4).*yy)) + exp(2*pi*i*F*(cos(alpha+2*pi/4).*xx+sin(alpha+2*pi/4).*yy))+ exp(2*pi*i*F*(cos(alpha+3*pi/4).*xx+sin(alpha+3*pi/4).*yy));
            phi3 = (real(phi3)/4);
            phi(1:end/2,1+end/2:end) = phi3(1:end/2,1+end/2:end);
            pos = find(phi<0);
            phi(pos) = 0;
            phi = phi/max(max(abs(phi)));
            phi = phi - mean(phi(:));
            phi(end/2+1:end,end/2+1:end) = -0.2;
            
        case 2
            sigma = 0.2; % parameter for spectral clustering
            deformFactor = 0.25;
            distType = 4;
            numAgl = 12;
            clsThre = 0.5;
            isScale = 1; %or 0, both work
            
            N = 512;
            x=0:1/N:(N-1)/N;
            [xo yo] = ndgrid(x);
            amp = 0.0;
            xx = xo + amp*sin(2*pi*xo);
            yy = yo + amp*sin(2*pi*yo);
            bmp = 0.0;
            alpha = -pi/8 + bmp*sin(2*pi*xo) + bmp*cos(2*pi*yo);
            F = 50;
            phi = exp(2*pi*i*F*(cos(alpha).*xx+sin(alpha).*yy)) + exp(2*pi*i*F*(cos(alpha+pi/3).*xx+sin(alpha+pi/3).*yy)) + exp(2*pi*i*F*(cos(alpha+2*pi/3).*xx+sin(alpha+2*pi/3).*yy));
            phi = real(phi)/3;
            % pause;
            alpha = -pi/8  + bmp*cos(2*pi*xo) + bmp*sin(2*pi*yo);
            F = 50;
            phi2 = exp(2*pi*i*F*(cos(alpha).*xx+sin(alpha).*yy)) + exp(2*pi*i*F*(cos(alpha+pi/2).*xx+sin(alpha+pi/2).*yy));
            phi2 = (real(phi2)/2);
            phi(end/2+1:end,1:end/2) = phi2(end/2+1:end,1:end/2);
            F = 60;
            phi3 = exp(2*pi*i*F*(cos(alpha).*xx+sin(alpha).*yy)) + exp(2*pi*i*F*(cos(alpha+pi/4).*xx+sin(alpha+pi/4).*yy)) + exp(2*pi*i*F*(cos(alpha+2*pi/4).*xx+sin(alpha+2*pi/4).*yy))+ exp(2*pi*i*F*(cos(alpha+3*pi/4).*xx+sin(alpha+3*pi/4).*yy));
            phi3 = (real(phi3)/4);
            phi(1:end/2,1+end/2:end) = phi3(1:end/2,1+end/2:end);
            pos = find(phi<0);
            phi(pos) = 0;
            phi = phi/max(max(abs(phi)));
            phi = phi - mean(phi(:));
            phi(end/2+1:end,end/2+1:end) = -0.2;
            
            alpha = pi/8 + bmp*sin(2*pi*xo) + bmp*cos(2*pi*yo);
            F = 90;
            phi4 = exp(2*pi*i*F*(cos(alpha).*xx+sin(alpha).*yy)) + exp(2*pi*i*F*(cos(alpha+pi/3).*xx+sin(alpha+pi/3).*yy)) + exp(2*pi*i*F*(cos(alpha+2*pi/3).*xx+sin(alpha+2*pi/3).*yy));
            phi4 = real(phi4)/3;
            % pause;
            alpha = pi/8  + bmp*cos(2*pi*xo) + bmp*sin(2*pi*yo);
            F = 90;
            phi2 = exp(2*pi*i*F*(cos(alpha).*xx+sin(alpha).*yy)) + exp(2*pi*i*F*(cos(alpha+pi/2).*xx+sin(alpha+pi/2).*yy));
            phi2 = (real(phi2)/2);
            phi4(end/2+1:end,1:end/2) = phi2(end/2+1:end,1:end/2);
            F = 100;
            phi3 = exp(2*pi*i*F*(cos(alpha).*xx+sin(alpha).*yy)) + exp(2*pi*i*F*(cos(alpha+pi/4).*xx+sin(alpha+pi/4).*yy)) + exp(2*pi*i*F*(cos(alpha+2*pi/4).*xx+sin(alpha+2*pi/4).*yy))+ exp(2*pi*i*F*(cos(alpha+3*pi/4).*xx+sin(alpha+3*pi/4).*yy));
            phi3 = (real(phi3)/4);
            phi4(1:end/2,1+end/2:end) = phi3(1:end/2,1+end/2:end);
            pos = find(phi4<0);
            phi4(pos) = 0;
            phi4 = phi4/max(max(abs(phi4)));
            phi4 = phi4 - mean(phi4(:));
            phi4(end/2+1:end,end/2+1:end) = -0.2;
            
            phi(end/4+1:3*end/4,end/4+1:3*end/4) = rot90(phi4(end/4+1:3*end/4,end/4+1:3*end/4),3);
            
            
        case 3
            sigma = 0.2; % parameter for spectral clustering
            deformFactor = 0.5;
            distType = 4;
            numAgl = 12;
            clsThre = 0.5;
            isScale = 1;
            
            N = 512;
            x=0:1/N:(N-1)/N;
            [xo yo] = ndgrid(x);
            amp = 0.0;
            xx = xo + amp*sin(2*pi*xo);
            yy = yo + amp*sin(2*pi*yo);
            bmp = 0.03;
            alpha = -pi/8 + bmp*sin(2*pi*xo) + bmp*cos(2*pi*yo);
            F = 50;
            phi = exp(2*pi*i*F*(cos(alpha).*xx+sin(alpha).*yy)) + exp(2*pi*i*F*(cos(alpha+pi/3).*xx+sin(alpha+pi/3).*yy)) + exp(2*pi*i*F*(cos(alpha+2*pi/3).*xx+sin(alpha+2*pi/3).*yy));
            phi = real(phi)/3;
            % pause;
            alpha = -pi/8  + bmp*cos(2*pi*xo) + bmp*sin(2*pi*yo);
            F = 50;
            phi2 = exp(2*pi*i*F*(cos(alpha).*xx+sin(alpha).*yy)) + exp(2*pi*i*F*(cos(alpha+pi/2).*xx+sin(alpha+pi/2).*yy));
            phi2 = (real(phi2)/2);
            phi(end/2+1:end,1:end/2) = phi2(end/2+1:end,1:end/2);
            F = 60;
            phi3 = exp(2*pi*i*F*(cos(alpha).*xx+sin(alpha).*yy)) + exp(2*pi*i*F*(cos(alpha+pi/4).*xx+sin(alpha+pi/4).*yy)) + exp(2*pi*i*F*(cos(alpha+2*pi/4).*xx+sin(alpha+2*pi/4).*yy))+ exp(2*pi*i*F*(cos(alpha+3*pi/4).*xx+sin(alpha+3*pi/4).*yy));
            phi3 = (real(phi3)/4);
            phi(1:end/2,1+end/2:end) = phi3(1:end/2,1+end/2:end);
            pos = find(phi<0);
            phi(pos) = 0;
            phi = phi/max(max(abs(phi)));
            phi = phi - mean(phi(:));
            phi(end/2+1:end,end/2+1:end) = -0.2;
            
        case 4
            sigma = 0.2; % parameter for spectral clustering
            deformFactor = 0.25;
            distType = 4;
            numAgl = 9;
            clsThre = 0.5;
            isScale = 1;
            
            N = patchSize;
            x=0:1/N:(N-1)/N;
            [xo yo] = ndgrid(x);
            amp = 0.0;
            xx = xo + amp*sin(2*pi*xo);
            yy = yo + amp*sin(2*pi*yo);
            bmp = 0.0;
            alpha = pi/12 + bmp*sin(2*pi*xo) + bmp*cos(2*pi*yo);
            F = 40;
            phi = exp(2*pi*i*F*(cos(alpha).*xx+sin(alpha).*yy)) + exp(2*pi*i*F*(cos(alpha+pi/3).*xx+sin(alpha+pi/3).*yy)) + exp(2*pi*i*F*(cos(alpha+2*pi/3).*xx+sin(alpha+2*pi/3).*yy));
            phi = real(phi)/3;
            pos = find(phi<0);
            phi(pos) = 0;
            
            
            phi2 = circshift(phi,3,1);
            phi2 = circshift(phi2,4,2);
            phi2 = max(phi2,phi);
            
            alpha = -pi/24  + bmp*cos(2*pi*xo) + bmp*sin(2*pi*yo);
            phi3 = exp(2*pi*i*F*(cos(alpha).*xx+sin(alpha).*yy)) + exp(2*pi*i*F*(cos(alpha+pi/2).*xx+sin(alpha+pi/2).*yy));
            phi3 = (real(phi3)/2);
            pos = find(phi3<0);
            phi3(pos) = 0;
            
            phi4 = circshift(phi3,4,1);
            phi4 = circshift(phi4,4,2);
            phi4 = max(phi4,phi3);
            
            phi(end/2+1:end,1:end/2) = phi2(1:end/2,1:end/2);
            phi(1:end/2,end/2+1:end) = phi3(1:end/2,1:end/2);
            phi(end/2+1:end,end/2+1:end) = phi4(1:end/2,1:end/2);
            
            phi(1:end/4,1:end/4) = -1.5*phi(1:end/4,1:end/4);
            phi(1:end/4,3*end/4+1:end) = phi(1:end/4,3*end/4+1:end)*2;
            phi(3*end/4+1:end,1:end/4) = phi(3*end/4+1:end,1:end/4)/2;
            phi(3*end/4+1:end,3*end/4+1:end) = -0.5*phi(3*end/4+1:end,3*end/4+1:end);
            
        case 5
            sigma = 0.2; % parameter for spectral clustering
            deformFactor = 0.5;
            distType = 4;
            numAgl = 9;
            clsThre = 0.5;
            isScale = 1;
            
            N = patchSize;
            x=0:1/N:(N-1)/N;
            [xo yo] = ndgrid(x);
            amp = 0.0;
            xx = xo + amp*sin(2*pi*xo);
            yy = yo + amp*sin(2*pi*yo);
            bmp = 0.0;
            alpha = pi/12 + bmp*sin(2*pi*xo) + bmp*cos(2*pi*yo);
            F = 40;
            phi = exp(2*pi*i*F*(cos(alpha).*xx+sin(alpha).*yy)) + exp(2*pi*i*F*(cos(alpha+pi/3).*xx+sin(alpha+pi/3).*yy)) + exp(2*pi*i*F*(cos(alpha+2*pi/3).*xx+sin(alpha+2*pi/3).*yy));
            phi = real(phi)/3;
            pos = find(phi<0);
            phi(pos) = 0;
            
            
            alpha = -pi/24  + bmp*cos(2*pi*xo) + bmp*sin(2*pi*yo);
            phi3 = exp(2*pi*i*F*(cos(alpha).*xx+sin(alpha).*yy)) + exp(2*pi*i*F*(cos(alpha+pi/2).*xx+sin(alpha+pi/2).*yy));
            phi3 = (real(phi3)/2);
            pos = find(phi3<0);
            phi3(pos) = 0;
            
            phi(:,end/2+1:end) = phi3(:,1:end/2);
            
            phi((1:end/3),1:end/2) = -phi((1:end/3),1:end/2);
            phi((1:end/3),end/2+1:end) = phi((1:end/3),end/2+1:end)*2;
            phi(2*end/3+1:end,1:end/2) = phi(2*end/3+1:end,1:end/2)/2;
            phi(2*end/3+1:end,end/2+1:end) = -phi(2*end/3+1:end,end/2+1:end);
    end
    
    phi = phi/max(max(abs(phi)));
    phi = phi - mean(phi(:));
    
    pic = figure;imagesc(phi);axis image;set(gca,'xtick',[]);set(gca,'ytick',[]);colormap (gray);
    set(gca, 'FontSize', 16);
    b=get(gca);
    set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
    head = sprintf('./results/clsImg%d.fig',ex);
    saveas(pic,head);
    head = sprintf('./results/clsImg%d.png',ex);
    saveas(pic,head);
    
    %set up parameters for ss transform
    para = struct('sigma',sigma,'gdSzType',gdSzType,'div',div,'NB',NB,'N',patchSize,'red',red,'epsl',1e-4,'is_real',...
        0,'t_sc',t_sc,'s_sc',s_sc,'num_direction',1,'deformFactor',deformFactor,'numAgl',numAgl,'distType',distType,'clsThre',clsThre,...
        'R_low2',R_low2,'R_low3',R_low3,'R_high2',R_high2,'R_high3',R_high3,'SPg',[],'spExtention',spExtention,'extention',extention,...
        'rad',rad,'threBD',threBD,'energyThre',energyThre,'isLake',isLake,'isScale',isScale,'typePick',2,'ex',ex);
    
    numAtom = 25;
    fudgeFactor = .5;
    % Estimate basis information, generate filtered images
    tic;
    [~,info,para] = crystalInfo(phi,para,is_show,numAtom,extention,fudgeFactor,fqThre*2);
    toc
    [typeMat,refType,R_high,refPos,timeRep,timeCls,timeSSTt] = skeletonCls(phi,para);
    timeRep
    timeCls
    timeSSTt
    
%     % close all
%     pic = figure;imagesc(skeletonShow(typeMat{1},512));axis image;set(gca,'xtick',[]);set(gca,'ytick',[]);colorbar;colormap (1-gray);
%     set(gca, 'FontSize', 16);
%     b=get(gca);
%     set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
%     head = sprintf('./results/clsTypeFnl%d.fig',ex);
%     saveas(pic,head);
%     head = sprintf('./results/clsTypeFnl%d.png',ex);
%     saveas(pic,head);
    
    if (ex == 1) || (ex == 3) || (ex == 4) % to get the successful rate of other examples, please modify the reference locations below
        szM = size(typeMat{1});
        numErr = 0;
        for cnti = 1:2
            for cntj = 1:2
                typeMat{1}(round(3*end/8)+szM(1)*(cnti-1)/2,round(3*end/8)+szM(2)*(cntj-1)/2)
                pos = find(typeMat{1}((1:end/2)+szM(1)*(cnti-1)/2,(1:end/2)+szM(2)*(cntj-1)/2)~= typeMat{1}(round(3*end/8)+szM(1)*(cnti-1)/2,round(3*end/8)+szM(2)*(cntj-1)/2));
                numErr = numErr + numel(pos);
            end
        end
        sucessRate = 1-numErr/prod(szM)
    end
    close all;
end

