function results = mainAnalysisSST(phiwhole,paraSST,info,is_show)
% This is the main code for analyzing crystal images
%
% Input:
% phiwhole: crystal image to be analyzed

%% Set up parameters and divide the whole image into small patches to fit computer RAM
SPg = paraSST.SPg;
[m,n] = size(phiwhole);
% Increase the size to divide the image into several pieces
phi = zeros(ceil(m/(0.5*paraSST.N))*0.5*paraSST.N,ceil(n/(0.5*paraSST.N))*0.5*paraSST.N);
[mm nn] = size(phi);

range1 = (1:m)+floor((mm-m)/2);
range2 = (1:n)+floor((nn-n)/2);
phi(range1,range2) = phiwhole;
[rm rn] = size(phi);
phiwhole = phi;
[m,n] = size(phiwhole);
numm = m/0.5/paraSST.N-1;
numn = n/0.5/paraSST.N-1;

BDCell = cell(numm,numn);
oriCell = cell(numm,numn);
maskCell = cell(numm,numn);
typeCell = cell(numm,numn);
lakeCell = cell(numm,numn);
spaceCell = cell(numm,numn);
GCell = cell(numm,numn);

%% Main loop for the analysis of each patch
[numm numn]
for cntm = 1:numm
    for cntn = 1:numn
        [cntm cntn]
        phi = phiwhole((1:paraSST.N)+(cntm-1)*paraSST.N/2,(1:paraSST.N)+(cntn-1)*paraSST.N/2);
        % Zero padding for non-periodic image
        [m n] = size(phi);
        dm = round(m*0.1/paraSST.div)*paraSST.div;
        dn = round(n*0.1/paraSST.div)*paraSST.div;
        phi = [zeros(dm,dn*2+n);zeros(m,dn) phi zeros(m,dn);zeros(dm,dn*2+n)];
        phi = -phi;
        
        % Clean peaks in the frequency domain
        [mphi nphi] = size(phi);
        [xg yg] = ndgrid(-mphi/2:mphi/2-1,-nphi/2:nphi/2-1);
        R = sqrt(xg.^2+yg.^2);
        Y = fftshift(fft2(phi));
        pos = find(R>paraSST.R_high3);
        Y(pos) = 0;
        pos = find(R<paraSST.R_low3);
        Y(pos) = 0;
        phi = real(ifft2(ifftshift(Y)));
        
        switch paraSST.isCheckType
            % 0: choose from cubic or hex lattice
            % 2: cubic
            % 3: hex
            case 0
                % Search for peaks
                num_wave = 3;%determined by the hexegonal reference grain
                R_low = paraSST.R_low3*(1-paraSST.spExtention); R_high = paraSST.R_high3*(1+paraSST.spExtention);
                %ss transform
                [ss_energy,ss_avgdx,ss_avgdy] = SS_ct2_polar_v2(paraSST.num_direction,phi,SPg,paraSST.NB,paraSST.rad,paraSST.is_real,[R_low,R_high],[0,2*pi],paraSST.epsl,paraSST.red,paraSST.t_sc,paraSST.s_sc);
                [~,TTEng_1st,TTEng_2nd] = LocSmooth(ss_energy,num_wave); % bin
                [agl,radii,~,~] = LocWeight(ss_energy,ss_avgdx,ss_avgdy,num_wave); % weight average
                %TTEng_1st = sum(TTEng_1st,3);
                %TTEng_2nd = sum(TTEng_2nd,3);
                
                num_wave = 2;%determined by the cubic reference grain
                R_low =paraSST.R_low2*(1-paraSST.spExtention); R_high = paraSST.R_high2*(1+paraSST.spExtention);
                %ss transform
                if (paraSST.R_low2~=paraSST.R_low3 | paraSST.R_high2~=paraSST.R_high3)
                    [ss_energy,ss_avgdx,ss_avgdy] = SS_ct2_polar_v2(paraSST.num_direction,phi,SPg,paraSST.NB,paraSST.rad,paraSST.is_real,[R_low,R_high],[0,2*pi],paraSST.epsl,paraSST.red,paraSST.t_sc,paraSST.s_sc);
                end
                [~,TTEng_1st2,TTEng_2nd2] = LocSmooth(ss_energy,num_wave);
                [agl2,radii2,~,~] = LocWeight(ss_energy,ss_avgdx,ss_avgdy,num_wave);
                %clear ss_energy ss_avgdx ss_avgdy
                %TTEng_1st2 = sum(TTEng_1st2,3);
                %TTEng_2nd2 = sum(TTEng_2nd2,3);
                
                % Determine the type of Bravais lattice
                posType = find(smoothImage(TTEng_1st2./(TTEng_1st2+TTEng_2nd2),5,1)>smoothImage(TTEng_1st./(TTEng_1st+TTEng_2nd),5,1));
                %posType = find(smoothImage(TTEng_1st2,5,1)>smoothImage(TTEng_1st,5,1));
                [szm,szn] = size(TTEng_1st);
                matType = 3*ones(szm,szn);
                matType(posType) = 2;
                
                % Correct small piece of cubic lattice
                binaryIndex = zeros(szm,szn);
                binaryIndex(posType) = 1;
                [CC,~] = periodicConnComp( 1-binaryIndex, 0 );
                %info.numPix = 100; check this
                for cnt = 1:CC.NumObjects
                    if numel(CC.PixelIdxList{cnt}) < info.numPix/(mphi/szm)/(nphi/szn)
                        binaryIndex(CC.PixelIdxList{cnt}) = 0;
                    end
                end
                posType3 = find(binaryIndex==0);
                matType(posType3) = 3;
                
                % Correct small piece of hexogonal lattice
                binaryIndex = zeros(szm,szn);
                posType3 = find(matType==3);
                binaryIndex(posType3) = 1;
                [CC,~] = periodicConnComp( 1-binaryIndex, 0 );
                for cnt = 1:CC.NumObjects
                    if numel(CC.PixelIdxList{cnt}) < info.numPix/(mphi/szm)/(nphi/szn)
                        binaryIndex(CC.PixelIdxList{cnt}) = 0;
                    end
                end
                posType2 = find(binaryIndex==0);
                matType(posType2) = 2;
                
                if paraSST.isCompDeform
                    % assumes y-coordinate pointing into direction of growing first index
                    agl = mod(agl+pi/18,pi)-pi/18;
                    pos = find(matType==3);
                    temp1 = radii(:,:,1);
                    temp2 = radii(:,:,2);
                    temp3 = radii(:,:,3);
                    avgRadii = median([temp1(pos); temp2(pos); temp3(pos)]);
                    % bring wave vectors into order according to angle
                    %                 [agl,perm] = sort(agl,3);
                    %                 for k = 1:size(perm,1)
                    %                     for l = 1:size(perm,2)
                    %                         radii(k,l,:) = radii(k,l,perm(k,l,:));
                    %                     end
                    %                 end
                    G = initGFromSSPeaks( agl, radii, pi/3*[0;1;2], avgRadii*ones(3,1));
                    %clear agl radii
                    
                    % assumes y-coordinate pointing into direction of growing first index
                    agl2 = mod(agl2+pi/18,pi)-pi/18;
                    if ~paraSST.isRectangle
                        % bring wave vectors into order according to angle
                        %                     [agl2,perm] = sort(agl2,3);
                        %                     for k = 1:size(perm,1)
                        %                         for l = 1:size(perm,2)
                        %                             radii2(k,l,:) = radii2(k,l,perm(k,l,:));
                        %                         end
                        %                     end
                        pos = find(matType==2);
                        temp1 = radii2(:,:,1);
                        temp2 = radii2(:,:,2);
                        avgRadii = median( [temp1(pos); temp2(pos)] );
                        G2 = initGFromSSPeaks( agl2, radii2, pi/2*[0;1], avgRadii*ones(2,1));
                    else
                        % bring wave vectors into order according to radius
                        [radii2,perm] = sort(radii2,3);
                        for k = 1:size(perm,1)
                            for l = 1:size(perm,2)
                                agl2(k,l,:) = agl2(k,l,perm(k,l,:));
                            end
                        end
                        temp1 = radii2(:,:,1);
                        avgRadii1 = median(temp1(:));
                        temp1 = radii2(:,:,2);
                        avgRadii2 = median(temp1(:));
                        G = initGFromSSPeaks( agl2, radii2, pi/2*[0;1], [avgRadii1;avgRadii2].*ones(2,1));
                    end
                    %clear avgRadii temp1 temp2 temp3 agl2 radii2
                    
                    % Choose deformation, SS energy, and angles
                    [pos1,pos2] = find(matType==2);
                    for cnt = 1:length(pos1)
                        G(pos1,pos2,:) = G2(pos1,pos2,:);
                    end
                    %clear G2
                    posType = find(matType == 2);
                    TTEng_1st(posType) = TTEng_1st2(posType);
                    TTEng_2nd(posType) = TTEng_2nd2(posType);
                    ori = zeros(szm,szn);
                    for cnt = 1:szm
                        for cnt2 = 1:szn
                            temp = [G(cnt,cnt2,1,1) G(cnt,cnt2,1,2);G(cnt,cnt2,2,1) G(cnt,cnt2,2,2)];
                            [U,H] = polarDec(temp);
                            if matType(cnt,cnt2)==3
                                ori(cnt,cnt2) = mod(real(acos(U(1,1))),pi/3)*180/pi;
                            else
                                ori(cnt,cnt2) = mod(pi-real(acos(U(1,1))),pi/2)*180/pi;
                            end
                        end
                    end
                else
                    for cnt = 1:szm
                        for cnt2 = 1:szn
                            if matType(cnt,cnt2)==3
                                ori(cnt,cnt2) = mod(median(agl(cnt,cnt2,:)),pi/3)*180/pi;
                            else
                                ori(cnt,cnt2) = mod(median(agl(cnt,cnt2,:)),pi/2)*180/pi;
                            end
                        end
                    end
                end
            case 2
                num_wave = 2;%determined by the cubic reference grain
                R_low =paraSST.R_low2*(1-paraSST.spExtention); R_high = paraSST.R_high2*(1+paraSST.spExtention);
                %ss transform
                [ss_energy,ss_avgdx,ss_avgdy] = SS_ct2_polar_v2(paraSST.num_direction,phi,SPg,paraSST.NB,paraSST.rad,paraSST.is_real,[R_low,R_high],[0,2*pi],paraSST.epsl,paraSST.red,paraSST.t_sc,paraSST.s_sc);
                [~,TTEng_1st,TTEng_2nd] = LocSmooth(ss_energy,num_wave);
                [agl,radii,~,~] = LocWeight(ss_energy,ss_avgdx,ss_avgdy,num_wave);
                %clear ss_energy ss_avgdx ss_avgdy
                
                % Determine the type of Bravais lattice
                [szm,szn] = size(TTEng_1st);
                matType = 2*ones(szm,szn);
                
                if paraSST.isCompDeform
                    % assumes y-coordinate pointing into direction of growing first index
                    agl = mod(agl+pi/18,pi)-pi/18;
                    if ~paraSST.isRectangle
                        % bring wave vectors into order according to angle
                        %                     [agl,perm] = sort(agl,3);
                        %                     for k = 1:size(perm,1)
                        %                         for l = 1:size(perm,2)
                        %                             radii(k,l,:) = radii(k,l,perm(k,l,:));
                        %                         end
                        %                     end
                        avgRadii = median(radii(:));
                        G = initGFromSSPeaks( agl, radii, pi/2*[0;1], avgRadii*ones(2,1));
                    else
                        % bring wave vectors into order according to radius
                        [radii,perm] = sort(radii,3);
                        for k = 1:size(perm,1)
                            for l = 1:size(perm,2)
                                agl(k,l,:) = agl(k,l,perm(k,l,:));
                            end
                        end
                        temp1 = radii(:,:,1);
                        avgRadii1 = median(temp1(:));
                        temp1 = radii(:,:,2);
                        avgRadii2 = median(temp1(:));
                        G = initGFromSSPeaks( agl, radii, pi/2*[0;1], [avgRadii1;avgRadii2].*ones(2,1));
                    end
                    %clear avgRadii temp1 temp2 temp3 agl2 radii2
                    
                    ori = zeros(szm,szn);
                    for cnt = 1:szm
                        for cnt2 = 1:szn
                            temp = [G(cnt,cnt2,1,1) G(cnt,cnt2,1,2);G(cnt,cnt2,2,1) G(cnt,cnt2,2,2)];
                            [U,H] = polarDec(temp);
                            ori(cnt,cnt2) = mod(pi-real(acos(U(1,1))),pi/2)*180/pi;
                        end
                    end
                else
                    for cnt = 1:szm
                        for cnt2 = 1:szn
                            ori(cnt,cnt2) = mod(median(agl(cnt,cnt2,:)),pi/2)*180/pi;
                        end
                    end
                end
            case 3
                % Search for peaks
                num_wave = 3;%determined by the hexegonal reference grain
                R_low = paraSST.R_low3*(1-paraSST.spExtention); R_high = paraSST.R_high3*(1+paraSST.spExtention);
                %ss transform
                [ss_energy,ss_avgdx,ss_avgdy] = SS_ct2_polar_v2(paraSST.num_direction,phi,SPg,paraSST.NB,paraSST.rad,paraSST.is_real,[R_low,R_high],[0,2*pi],paraSST.epsl,paraSST.red,paraSST.t_sc,paraSST.s_sc);
                [~,TTEng_1st,TTEng_2nd] = LocSmooth(ss_energy,num_wave); % bin
                [agl,radii,~,~] = LocWeight(ss_energy,ss_avgdx,ss_avgdy,num_wave); % weight average
                %TTEng_1st = sum(TTEng_1st,3);
                %TTEng_2nd = sum(TTEng_2nd,3);
                
                % Determine the type of Bravais lattice
                [szm,szn] = size(TTEng_1st);
                matType = 3*ones(szm,szn);
                
                if paraSST.isCompDeform
                    % assumes y-coordinate pointing into direction of growing first index
                    agl = mod(agl+pi/18,pi)-pi/18;
                    pos = find(matType==3);
                    temp1 = radii(:,:,1);
                    temp2 = radii(:,:,2);
                    temp3 = radii(:,:,3);
                    avgRadii = median([temp1(pos); temp2(pos); temp3(pos)]);
                    % bring wave vectors into order according to angle
                    %                 [agl,perm] = sort(agl,3);
                    %                 for k = 1:size(perm,1)
                    %                     for l = 1:size(perm,2)
                    %                         radii(k,l,:) = radii(k,l,perm(k,l,:));
                    %                     end
                    %                 end
                    G = initGFromSSPeaks( agl, radii, pi/3*[0;1;2], avgRadii*ones(3,1));
                    %clear agl radii
                    
                    ori = zeros(szm,szn);
                    for cnt = 1:szm
                        for cnt2 = 1:szn
                            temp = [G(cnt,cnt2,1,1) G(cnt,cnt2,1,2);G(cnt,cnt2,2,1) G(cnt,cnt2,2,2)];
                            [U,H] = polarDec(temp);
                            ori(cnt,cnt2) = mod(real(acos(U(1,1))),pi/3)*180/pi;
                        end
                    end
                else
                    for cnt = 1:szm
                        for cnt2 = 1:szn
                            ori(cnt,cnt2) = mod(median(agl(cnt,cnt2,:)),pi/3)*180/pi;
                        end
                    end
                end
        end
        
        % Compute grain boundary
        switch paraSST.typeBD
            case 2
                BD = 1./sqrt(2-TTEng_2nd./TTEng_1st);
            case 1
                BD = 1./sqrt(TTEng_1st-TTEng_2nd+1);
        end
        valBD = max(BD(:));
        
        % Enhance the intensity of isolated defects
        % Need to do it locally, because this is not robust to
        % global change.
        % Cannot distinguish full grain and no grain
        posBD = find(TTEng_1st+TTEng_2nd>max(TTEng_1st(:)+TTEng_2nd(:))/10);
        vec = TTEng_1st(posBD)+TTEng_2nd(posBD);
        posBD = find(TTEng_1st+TTEng_2nd<median(vec(:))*paraSST.threBD);
        BD(posBD) = valBD;
        
        % Compute a mask area indicating uncertain area including
        % isolated defects in this step
        mask = zeros(size(BD));
        mask(posBD) = 1;
        
        % Cut the zero padding around the non-periodic image
        BD = BD(dm/paraSST.div+1:dm/paraSST.div+m/paraSST.div,dn/paraSST.div+1:dn/paraSST.div+n/paraSST.div);
        BDCell{cntm,cntn} = BD;
        
        ori = ori(dm/paraSST.div+1:dm/paraSST.div+m/paraSST.div,dn/paraSST.div+1:dn/paraSST.div+n/paraSST.div);
        oriCell{cntm,cntn} = ori;
        
        mask = mask(dm/paraSST.div+1:dm/paraSST.div+m/paraSST.div,dn/paraSST.div+1:dn/paraSST.div+n/paraSST.div);
        
        matType = matType(dm/paraSST.div+1:dm/paraSST.div+m/paraSST.div,dn/paraSST.div+1:dn/paraSST.div+n/paraSST.div);
        
        if paraSST.isCompDeform
            G1 = zeros([size(mask),2,2]);
            for cnt1 = 1:2
                for cnt2 = 1:2
                    G1(:,:,cnt1,cnt2) = G(dm/paraSST.div+1:dm/paraSST.div+m/paraSST.div,dn/paraSST.div+1:dn/paraSST.div+n/paraSST.div,cnt1,cnt2);
                end
            end
        end
        
        % Correct small piece of mask, excluding small isolated defetcs,
        % only contains the liquid part
        [CC,~] = periodicConnComp(1-mask,0);
        for cnt = 1:CC.NumObjects
            if numel(CC.PixelIdxList{cnt}) < info.numPix/(mphi/szm)/(nphi/szn)%*16
                mask(CC.PixelIdxList{cnt}) = 0;
            end
        end
        posBD = find(mask>0);
        
        % Update the uncertain area
        %matType(posBD) = 1; % no need to do this. 12/29/2016
        maskCell{cntm,cntn} = mask;
        typeCell{cntm,cntn} = matType;
        if paraSST.isCompDeform
            GCell{cntm,cntn} = G1;
        end
        %         if is_show
        %             pic = figure;imagesc(phi(dm+1:end-dm,dn+1:end-dn));axis square;title('Original image');set(gca,'xtick',[],'ytick',[]); colormap (1-gray);
        %             pic = figure;imagesc(BD);axis square;title('Detected defect boundary of type IV');set(gca,'xtick',[],'ytick',[]); colormap (1-gray);
        %             pic = figure;imagesc(ori); axis square;title('Orientation via weighted mean');set(gca,'xtick',[],'ytick',[]);caxis([0, 60]);colorbar;colormap(pmkmp());
        %             pic = figure;imagesc(matType);colorbar;title('Type of bravais lattices');axis square;
        %             pause;
        %         end
        %         close all;
    end
end
%clear G G1

%% Assemble results from small patches
SPg(1) = length(dm/paraSST.div+1:dm/paraSST.div+m/paraSST.div);
SPg(2) = length(dn/paraSST.div+1:dn/paraSST.div+n/paraSST.div);
BD = zeros(SPg/2.*[numm+1 numn+1]);
ori = zeros(SPg/2.*[numm+1 numn+1]);
mask = zeros(SPg/2.*[numm+1 numn+1]);
matType = zeros(SPg/2.*[numm+1 numn+1]);
if paraSST.isCompDeform
    deformG= zeros([SPg/2.*[numm+1 numn+1],2,2]);
end
for cntm = 1:numm
    for cntn = 1:numn
        if cntm == 1 & cntn == 1
            temp = BDCell{cntm,cntn};
            BD(1:SPg(1),1:SPg(2)) = temp(1:SPg(1),1:SPg(2));
            temp = oriCell{cntm,cntn};
            ori(1:SPg(1),1:SPg(2)) = temp(1:SPg(1),1:SPg(2));
            temp = maskCell{cntm,cntn};
            mask(1:SPg(1),1:SPg(2)) = temp(1:SPg(1),1:SPg(2));
            temp = typeCell{cntm,cntn};
            matType(1:SPg(1),1:SPg(2)) = temp(1:SPg(1),1:SPg(2));
            if paraSST.isCompDeform
                temp = GCell{cntm,cntn};
                for cnt1 = 1:2
                    for cnt2 = 1:2
                        deformG(1:SPg(1),1:SPg(2),cnt1,cnt2) = temp(1:SPg(1),1:SPg(2),cnt1,cnt2);
                    end
                end
            end
        else if cntm ==1 & cntn ~= 1
                temp = BDCell{cntm,cntn};
                BD(1:SPg(1),((SPg(2)*0.25+1):SPg(2))+0.5*SPg(2)*(cntn-1)) = temp(1:SPg(1),SPg(2)*0.25+1:end);
                temp = oriCell{cntm,cntn};
                ori(1:SPg(1),((SPg(2)*0.25+1):SPg(2))+0.5*SPg(2)*(cntn-1)) = temp(1:SPg(1),SPg(2)*0.25+1:end);
                temp = maskCell{cntm,cntn};
                mask(1:SPg(1),((SPg(2)*0.25+1):SPg(2))+0.5*SPg(2)*(cntn-1)) = temp(1:SPg(1),SPg(2)*0.25+1:end);
                temp = typeCell{cntm,cntn};
                matType(1:SPg(1),((SPg(2)*0.25+1):SPg(2))+0.5*SPg(2)*(cntn-1)) = temp(1:SPg(1),SPg(2)*0.25+1:end);
                if paraSST.isCompDeform
                    temp = GCell{cntm,cntn};
                    for cnt1 = 1:2
                        for cnt2 = 1:2
                            deformG(1:SPg(1),((SPg(2)*0.25+1):SPg(2))+0.5*SPg(2)*(cntn-1),cnt1,cnt2) = temp(1:SPg(1),SPg(2)*0.25+1:end,cnt1,cnt2);
                        end
                    end
                end
            else if cntm ~= 1 & cntn == 1
                    temp = BDCell{cntm,cntn};
                    BD(((SPg(1)*0.25+1):SPg(1))+0.5*SPg(1)*(cntm-1),1:SPg(2)) = temp(SPg(1)*0.25+1:end,1:SPg(2));
                    temp = oriCell{cntm,cntn};
                    ori(((SPg(1)*0.25+1):SPg(1))+0.5*SPg(1)*(cntm-1),1:SPg(2)) = temp(SPg(1)*0.25+1:end,1:SPg(2));
                    temp = maskCell{cntm,cntn};
                    mask(((SPg(1)*0.25+1):SPg(1))+0.5*SPg(1)*(cntm-1),1:SPg(2)) = temp(SPg(1)*0.25+1:end,1:SPg(2));
                    temp = typeCell{cntm,cntn};
                    matType(((SPg(1)*0.25+1):SPg(1))+0.5*SPg(1)*(cntm-1),1:SPg(2)) = temp(SPg(1)*0.25+1:end,1:SPg(2));
                    if paraSST.isCompDeform
                        temp = GCell{cntm,cntn};
                        for cnt1 = 1:2
                            for cnt2 = 1:2
                                deformG(((SPg(1)*0.25+1):SPg(1))+0.5*SPg(1)*(cntm-1),1:SPg(2),cnt1,cnt2) = temp(SPg(1)*0.25+1:end,1:SPg(2),cnt1,cnt2);
                            end
                        end
                    end
                else
                    temp = BDCell{cntm,cntn};
                    BD(((SPg(1)*0.25+1):SPg(1))+0.5*SPg(1)*(cntm-1),((SPg(2)*0.25+1):SPg(2))+0.5*SPg(2)*(cntn-1)) = temp(SPg(1)*0.25+1:end,SPg(2)*0.25+1:end);
                    temp = oriCell{cntm,cntn};
                    ori(((SPg(1)*0.25+1):SPg(1))+0.5*SPg(1)*(cntm-1),((SPg(2)*0.25+1):SPg(2))+0.5*SPg(2)*(cntn-1)) = temp(SPg(1)*0.25+1:end,SPg(2)*0.25+1:end);
                    temp = maskCell{cntm,cntn};
                    mask(((SPg(1)*0.25+1):SPg(1))+0.5*SPg(1)*(cntm-1),((SPg(2)*0.25+1):SPg(2))+0.5*SPg(2)*(cntn-1)) = temp(SPg(1)*0.25+1:end,SPg(2)*0.25+1:end);
                    temp = typeCell{cntm,cntn};
                    matType(((SPg(1)*0.25+1):SPg(1))+0.5*SPg(1)*(cntm-1),((SPg(2)*0.25+1):SPg(2))+0.5*SPg(2)*(cntn-1)) = temp(SPg(1)*0.25+1:end,SPg(2)*0.25+1:end);
                    if paraSST.isCompDeform
                        temp = GCell{cntm,cntn};
                        for cnt1 = 1:2
                            for cnt2 = 1:2
                                deformG(((SPg(1)*0.25+1):SPg(1))+0.5*SPg(1)*(cntm-1),((SPg(2)*0.25+1):SPg(2))+0.5*SPg(2)*(cntn-1),cnt1,cnt2) = temp(SPg(1)*0.25+1:end,SPg(2)*0.25+1:end,cnt1,cnt2);
                            end
                        end
                    end
                end
            end
        end
    end
end
[m n] = size(BD);
BD = BD(  ceil(range1(1)*m/rm):floor(range1(end)*m/rm) , ceil(range2(1)*n/rn):floor(range2(end)*n/rn) );
ori = ori(  ceil(range1(1)*m/rm):floor(range1(end)*m/rm) , ceil(range2(1)*n/rn):floor(range2(end)*n/rn) );
mask = mask(ceil(range1(1)*m/rm):floor(range1(end)*m/rm) , ceil(range2(1)*n/rn):floor(range2(end)*n/rn) );
matType = matType(ceil(range1(1)*m/rm):floor(range1(end)*m/rm) , ceil(range2(1)*n/rn):floor(range2(end)*n/rn) );
if paraSST.isCompDeform
    tempG = deformG;
    deformG = zeros([size(BD),2,2]);
    for cnt1 = 1:2
        for cnt2 = 1:2
            deformG(:,:,cnt1,cnt2) = tempG(ceil(range1(1)*m/rm):floor(range1(end)*m/rm) , ceil(range2(1)*n/rn):floor(range2(end)*n/rn),cnt1,cnt2);
        end
    end
end
%clear tempG;
val = min(BD(:));
BD = BD-val;
val = max(BD(:));
%if val > 0.5
BD = BD/val;
%end

% if is_show
%     figure;imagesc(info.dphi1);colorbar;title('dphi1');caxis([0 1]);
%     figure;imagesc(info.dphi2);colorbar;title('dphi2');caxis([0 1]);
%     figure;imagesc(info.matLake);colorbar;title('matLake');caxis([0 1]);
% end

%% Image segmentation for lakes and outliers
% Combine the results by energy and variation
pos = find(matType==1);
temp = ones(size(matType));
temp(pos) = 0;
[szLake1 szLake2] = size(info.matLake);
temp = imresize(temp,[szLake1 szLake2]);
info.matLake = max(info.matLake,temp);% use max to correct the errors of each method
pos = find(info.matLake>0);
info.matLake = zeros([szLake1,szLake2]);
info.matLake(pos) = 1;
% if is_show
%     figure;imagesc(info.matLake);colorbar;title('matLake2');caxis([0 1]);
% end
info.matLake = min(info.dphi,info.matLake); % use min to combine results

outlier = zeros(szLake1,szLake2); % isolated atoms in lake or isolated small vacancy

% Correct small pieces of lake outliers in grain
info.matLake = smoothImage(info.matLake,10,3);
pos = find(info.matLake>0.1);
info.matLake = zeros(size(info.matLake));
info.matLake(pos) = 1;
if ~any(strcmp('numAtomOutlier',fieldnames(paraSST)))
    numAtomOutlier = 16;
else
    numAtomOutlier = paraSST.numAtomOutlier;
end
[CC,~] = periodicConnComp(info.matLake, 0 );
for cnt = 1:CC.NumObjects
    if numel(CC.PixelIdxList{cnt}) < numAtomOutlier*info.numPixAtom^2 % number of pixels for numAtomOutlier atoms
        info.matLake(CC.PixelIdxList{cnt}) = 1;
        outlier(CC.PixelIdxList{cnt}) = 1;
    end
end
% Before the following step, 0 for lake and 1 for grain
info.matLake = 1-info.matLake;
% After the above step, 0 for grain and 1 for lake
% Correct small pieces of grain outliers in lake
info.matLake = smoothImage(info.matLake,10,3);
pos = find(info.matLake>0.1);
info.matLake = zeros(size(info.matLake));
info.matLake(pos) = 1;


[CC,~] = periodicConnComp(info.matLake, 0 );
for cnt = 1:CC.NumObjects
    if numel(CC.PixelIdxList{cnt}) < numAtomOutlier*info.numPixAtom^2 % number of pixels for numAtomOutlier atoms
        info.matLake(CC.PixelIdxList{cnt}) = 1;
        outlier(CC.PixelIdxList{cnt}) = 1;
    end
end
[szBD1 szBD2] = size(BD);
info.matLake = imresize(info.matLake,[szBD1 szBD2]);
pos = find(info.matLake>1/10);
info.matLake = zeros(size(info.matLake));
info.matLake(pos) = 1;

% Update BD with info.matLake
if paraSST.isLake
    pos = find(info.matLake>0);
    val = max(BD(:));
    BD(pos) = val;
end
%% Image segmentation for grains
temp = edge(matType,'sobel');
posEdge = find(temp==1);
temp = BD;
pos= find(matType==2);
temp(pos) = 0;
temp(posEdge) = 1;
temp = smoothImage(3*(1-temp), 2, 1);
val = min(temp(:));
temp = temp-val;
val = max(temp(:));
temp = 3*temp/val;
binaryIndexBD = findDeftArea(temp,paraSST.energyThre,1,0);
pos = find(info.matLake>0);
binaryIndexBD(pos) = 1;
[CC,~] = periodicConnComp(binaryIndexBD, 0 );
numAtomOutlier = 16;
for cnt = 1:CC.NumObjects
    if numel(CC.PixelIdxList{cnt}) < numAtomOutlier*info.numPixAtom^2/(szLake1/szBD1)/(szLake2/szBD2)
        binaryIndexBD(CC.PixelIdxList{cnt}) = 1;
    end
end
[CC,~] = periodicConnComp(binaryIndexBD, 0 );
temp = zeros(szBD1,szBD2);
vec = round(rand(1,CC.NumObjects)*(10000-CC.NumObjects)) + (1:CC.NumObjects);
[vec,ord] = sort(vec);
for cnt = 1:CC.NumObjects
    temp(CC.PixelIdxList{cnt}) = ord(cnt);
end
segmentation = temp;

%% Summarize results
results = struct('BD',[],'ori',[],'SPg',[],'mask',[],'matType',[],'matLake',[],'outlier',[],'segmentation',[],'deformG',[]);
results.BD = BD;
results.ori = ori;
results.SPg = SPg;
results.mask = mask;
results.matType = matType;
results.matLake = info.matLake;
results.space = info.space;
results.outlier = outlier;
results.segmentation = segmentation;
if paraSST.isCompDeform
    results.deformG = deformG;
end
