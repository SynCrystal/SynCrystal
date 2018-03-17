N = 64*2;
x=0:1/N:(N-1)/N;
[xo yo] = ndgrid(x);

% generate data for fig 14
if 0
    amp = 0.0;
    
    for agl = 0:35
        xx = xo + amp*sin(2*pi*xo);
        yy = yo + amp*sin(2*pi*yo);
        bmp = 0.0;
        alpha = agl*pi/36 + bmp*sin(2*pi*xo) + bmp*cos(2*pi*yo);
        F = 10;
        phi = exp(2*pi*i*F*(cos(alpha).*xx+sin(alpha).*yy)) + exp(2*pi*i*F*(cos(alpha+pi/3).*xx+sin(alpha+pi/3).*yy)) + exp(2*pi*i*F*(cos(alpha+2*pi/3).*xx+sin(alpha+2*pi/3).*yy));
        phi = real(phi)/3;
        pos = find(phi<0);
        phi(pos) = 0;
        
        
        pic = figure;imagesc(phi(end-63:end,end-63:end)); axis square;axis off;
        str = sprintf('./data/synthetic/fig14/1/%d',agl);
        print(gcf, '-depsc2', str);
        
        
        phi2 = circshift(phi,3,1);
        phi2 = circshift(phi2,4,2);
        phi2 = max(phi2,phi);
        
        pic = figure;imagesc(phi2(end-63:end,end-63:end)); axis square;axis off;
        str = sprintf('./data/synthetic/fig14/2/%d',agl);
        print(gcf, '-depsc2', str);
        
        alpha = -pi/24  + bmp*cos(2*pi*xo) + bmp*sin(2*pi*yo);
        phi3 = exp(2*pi*i*F*(cos(alpha).*xx+sin(alpha).*yy)) + exp(2*pi*i*F*(cos(alpha+pi/2).*xx+sin(alpha+pi/2).*yy));
        phi3 = (real(phi3)/2);
        pos = find(phi3<0);
        phi3(pos) = 0;
        
        pic = figure;imagesc(phi3(end-63:end,end-63:end)); axis square;axis off;
        str = sprintf('./data/synthetic/fig14/3/%d',agl);
        print(gcf, '-depsc2', str);
        
        phi4 = circshift(phi3,4,1);
        phi4 = circshift(phi4,4,2);
        phi4 = max(phi4,phi3);
        
        pic = figure;imagesc(phi4(end-63:end,end-63:end)); axis square;axis off;
        str = sprintf('./data/synthetic/fig14/4/%d',agl);
        print(gcf, '-depsc2', str);
        
        close all;
    end
end

if 1
    N = 512;
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
    
    indm = ones(size(phi));
    indm(end/2+1:end,1:end/2) = 2;
    indm(1:end/2,end/2+1:end) = 3;
    indm(end/2+1:end,end/2+1:end) = 4;
    
    cnt1 = 0; cnt2 = 0; cnt3 = 0; cnt4 = 0;
    for i1 = 1:16:N
        for i2 = 1:16:N
            type = indm(i1,i2);
            rd = i1-31:i1+32; cl = i2-31:i2+32;
            img = phi(mod(rd-1,N)+1,mod(cl-1,N)+1);
            pic = figure;imagesc(img); axis square;axis off;
            switch type
                case 1
                    str = sprintf('./data/synthetic/fig14/1/%d',cnt1+36); cnt1 = cnt1 + 1;
                case 2
                    str = sprintf('./data/synthetic/fig14/2/%d',cnt2+36); cnt2 = cnt2 + 1;
                case 3
                    str = sprintf('./data/synthetic/fig14/3/%d',cnt3+36); cnt3 = cnt3 + 1;
                case 4
                    str = sprintf('./data/synthetic/fig14/4/%d',cnt4+36); cnt4 = cnt4 + 1;
            end
            print(gcf, '-depsc2', str);
            close all;
        end
    end
end

% generate data for fig 12
if 0
    amp = 0.0;
    xx = xo + amp*sin(2*pi*xo);
    yy = yo + amp*sin(2*pi*yo);
    
    for agl = 0:35
        bmp = 0.03;
        alpha = agl*pi/36 + bmp*sin(2*pi*xo) + bmp*cos(2*pi*yo);
        F = 12.5;
        phi = exp(2*pi*i*F*(cos(alpha).*xx+sin(alpha).*yy)) + exp(2*pi*i*F*(cos(alpha+pi/3).*xx+sin(alpha+pi/3).*yy)) + exp(2*pi*i*F*(cos(alpha+2*pi/3).*xx+sin(alpha+2*pi/3).*yy));
        phi = real(phi)/3;
        
        pic = figure;imagesc(phi(end-63:end,end-63:end)); axis square;axis off;
        str = sprintf('./data/synthetic/fig12/1/%d',agl);
        print(gcf, '-depsc2', str);
        
        alpha = agl*pi/36  + bmp*cos(2*pi*xo) + bmp*sin(2*pi*yo);
        F = 12.5;
        phi2 = exp(2*pi*i*F*(cos(alpha).*xx+sin(alpha).*yy)) + exp(2*pi*i*F*(cos(alpha+pi/2).*xx+sin(alpha+pi/2).*yy));
        phi2 = (real(phi2)/2);
        
        pic = figure;imagesc(phi2(end-63:end,end-63:end)); axis square;axis off;
        str = sprintf('./data/synthetic/fig12/2/%d',agl);
        print(gcf, '-depsc2', str);
        
        F = 15;
        phi3 = exp(2*pi*i*F*(cos(alpha).*xx+sin(alpha).*yy)) + exp(2*pi*i*F*(cos(alpha+pi/4).*xx+sin(alpha+pi/4).*yy)) + exp(2*pi*i*F*(cos(alpha+2*pi/4).*xx+sin(alpha+2*pi/4).*yy))+ exp(2*pi*i*F*(cos(alpha+3*pi/4).*xx+sin(alpha+3*pi/4).*yy));
        phi3 = (real(phi3)/4);
        
        pic = figure;imagesc(phi3(end-63:end,end-63:end)); axis square;axis off;
        str = sprintf('./data/synthetic/fig12/3/%d',agl);
        print(gcf, '-depsc2', str);
        
        pic = figure;imagesc(ones(64,64)*randn(1,1)*5); axis square;axis off;
        str = sprintf('./data/synthetic/fig12/4/%d',agl);
        print(gcf, '-depsc2', str);
        
        close all;
    end
end

if 1
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
    
    indm = ones(size(phi));
    indm(end/2+1:end,1:end/2) = 2;
    indm(1:end/2,1+end/2:end) = 3;
    indm(end/2+1:end,end/2+1:end) = 4;
    
    cnt1 = 0; cnt2 = 0; cnt3 = 0; cnt4 = 0;
    for i1 = 1:16:N
        for i2 = 1:16:N
            type = indm(i1,i2);
            rd = i1-31:i1+32; cl = i2-31:i2+32;
            img = phi(mod(rd-1,N)+1,mod(cl-1,N)+1);
            pic = figure;imagesc(img); axis square;axis off;
            switch type
                case 1
                    str = sprintf('./data/synthetic/fig12/1/%d',cnt1+36); cnt1 = cnt1 + 1;
                case 2
                    str = sprintf('./data/synthetic/fig12/2/%d',cnt2+36); cnt2 = cnt2 + 1;
                case 3
                    str = sprintf('./data/synthetic/fig12/3/%d',cnt3+36); cnt3 = cnt3 + 1;
                case 4
                    str = sprintf('./data/synthetic/fig12/4/%d',cnt4+36); cnt4 = cnt4 + 1;
            end
            print(gcf, '-depsc2', str);
            close all;
        end
    end
end
