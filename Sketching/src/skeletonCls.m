function [typeMatCell,refSktCell,R_high,refPos,timeRep,timeCls,timeSSTt] = skeletonCls(phi,para)
% This is the main code for analyzing crystal images
%
% Input:
% phi: crystal image to be analyzed

%% Set up parameters and divide the whole image into small patches to fit computer RAM
sigma = para.sigma;
gdSzType = para.gdSzType;
[m,n] = size(phi);
% Increase the size to divide the image into several pieces
phiTemp = zeros(ceil(m/(0.5*para.N))*0.5*para.N,ceil(n/(0.5*para.N))*0.5*para.N);
[mm nn] = size(phiTemp);

range1 = (1:m)+floor((mm-m)/2);
range2 = (1:n)+floor((nn-n)/2);
phiTemp(range1,range2) = phi;
[rm rn] = size(phiTemp);
phi = phiTemp;
[m,n] = size(phi);
numm = m/0.5/para.N-1;
numn = n/0.5/para.N-1;

refSktCell = cell(numm,numn);
typeMatCell = cell(numm,numn);
timeRep = 0;
timeSSTt = 0;
timeCls = 0;
%% Main loop for the analysis of each patch
[numm numn]
for cntm = 1:numm
    for cntn = 1:numn
        %[cntm cntn]
        
        phiTemp = phi((1:para.N)+(cntm-1)*para.N/2,(1:para.N)+(cntn-1)*para.N/2);
        % Zero padding for non-periodic image
        [m n] = size(phiTemp);
        dm = 0;
        dn = 0;
        phiTemp = [zeros(dm,dn*2+n);zeros(m,dn) phiTemp zeros(m,dn);zeros(dm,dn*2+n)];
        phiTemp = -phiTemp;
        
        % Clean peaks in the frequency domain
        [mphi nphi] = size(phiTemp);
        [xg yg] = ndgrid(-mphi/2:mphi/2-1,-nphi/2:nphi/2-1);
        R = sqrt(xg.^2+yg.^2);
        Y = fftshift(fft2(phiTemp));
        pos = find(R>para.R_high3);
        Y(pos) = 0;
        pos = find(R<para.R_low3);
        Y(pos) = 0;
        phiTemp = real(ifft2(ifftshift(Y)));
        
        % Search for peaks
        R_low = min(para.R_low2,para.R_low3)*(1-para.spExtention); R_high = 2*max(para.R_high2,para.R_high3)*(1+para.spExtention);
        [skeletonWhole,~,~,~,~,~,~,timeSST,timeSkt]  = skeletonPolar(para.num_direction,phiTemp,para.SPg,para.NB,para.rad,para.is_real,...
            R_low,R_high,para.epsl,para.red,para.t_sc,para.s_sc,para.deformFactor,para.numAgl,para.isScale);
        szSSTWhole = size(skeletonWhole);
        pos1 = gdSzType:gdSzType:szSSTWhole(3)-gdSzType; pos2 = gdSzType:gdSzType:szSSTWhole(4)-gdSzType;
        skeleton = skeletonWhole(:,:,pos1,pos2);
        szSST = size(skeleton);
        timeRep = timeRep + timeSkt;
        timeSSTt = timeSSTt + timeSST;
        
        figure;imagesc(skeleton(:,:,floor(szSST(3)/2), floor(szSST(4)/2))); colormap(1-gray);axis xy; colorbar;
        
        tic;
        skeleton = reshape(skeleton,[szSST(1),szSST(2),szSST(3)*szSST(4)]);
        [adjMat,~] = AdjMatSkeleton(skeleton,sigma,para.distType,para.clsThre);
        [group_resid , ~, ~ , num_group_resid] = SpectralClustering_est(adjMat,3);
        typeMatCell{cntm,cntn} = zeros(szSST(3)*szSST(4),1);
        rec = [];
        tp = 0;
        refCell = cell(1,num_group_resid);
        for cnt = 1:num_group_resid
            posg = find(group_resid==cnt);
            % compute average skeleton
            refCell{tp+1} = median(skeleton(:,:,posg),3);
            if numel(posg) > sqrt(szSST(3)*szSST(4))
                tp = tp + 1;
                rec = [rec,posg(1)];
                typeMatCell{cntm,cntn}(posg) = tp;
                group_resid(posg) = tp;
            else
                group_resid(posg) = 0;
            end
        end
        
        % get the best reference skeleton
        group_resid = reshape(group_resid,[szSST(3),szSST(4)]);
        refPos = zeros(1,tp);
        for cnt = 1:tp
            temp = zeros(szSST(3),szSST(4));
            pos = find(group_resid==cnt);
            temp(pos) = 1;
            temp = smoothImage(temp, round(szSST(3)/4), szSST(3)/8);
            [~,pos] = max(temp(:));
            refCell{cnt} = skeleton(:,:,pos);
            refPos(cnt) = pos;
        end
        refSktCell{cntm,cntn} = cell(1,tp);
        for cnt = 1:tp
            refSktCell{cntm,cntn}{cnt} = refCell{cnt};
        end
        typeMatCell{cntm,cntn} = reshape(typeMatCell{cntm,cntn},[szSST(3),szSST(4)]);
        timeCls = timeCls + toc;
        if 1
            pic = figure;imagesc(skeletonShow(typeMatCell{cntm,cntn},512));axis image;set(gca,'xtick',[]);set(gca,'ytick',[]);colorbar;colormap (1-gray);
            hold on; sz = size(typeMatCell{cntm,cntn});
            for cnt = 1:numel(refPos)
                [idx1,idx2] = ind2sub(sz,refPos(cnt));
                plot(idx2/sz(2)*512,idx1/sz(1)*512,'r.','MarkerSize',40);
            end
            set(gca, 'FontSize', 16);
            b=get(gca);
            set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
            head = sprintf('./results/clsType%d.fig',para.ex);
            saveas(pic,head);
            head = sprintf('./results/clsType%d.png',para.ex);
            saveas(pic,head);
            
            for cnt = 1:length(refSktCell{cntm,cntn})
                pic = figure;imagesc([0,pi],[0,R_high/2],skeletonShow(refSktCell{cntm,cntn}{cnt}(1:end/2,:),512));axis square;xlabel('Angle');ylabel('Radius');colorbar;colormap (1-gray);caxis([0,1]); axis xy;
                set(gca, 'FontSize', 16);
                b=get(gca);
                set(b.XLabel, 'FontSize', 16);set(b.YLabel, 'FontSize', 16);set(b.ZLabel, 'FontSize', 16);set(b.Title, 'FontSize', 16);
                head = sprintf('./results/clsSSTType%d_%d.fig',para.ex,cnt);
                saveas(pic,head);
                head = sprintf('./results/clsSSTType%d_%d.png',para.ex,cnt);
                saveas(pic,head);
            end
        end
        
       % typeMatCell{cntm,cntn} =
       % skeletonPick(skeletonWhole,refSktCell{cntm,cntn},para); % refine
       % the result using identified reference lattice representation
    end
end