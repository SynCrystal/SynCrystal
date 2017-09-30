function [typeCell,typeMatCell,ss_energy,adjMat,ssEnergySpl] = crystalType(phiwhole,paraSST,is_show)
% This is the main code for analyzing crystal images
%
% Input:
% phiwhole: crystal image to be analyzed

%% Set up parameters and divide the whole image into small patches to fit computer RAM
SPg = paraSST.SPg;
sigma = paraSST.sigma;
gdSzType = paraSST.gdSzType;
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


typeCell = cell(numm,numn);
typeMatCell = cell(numm,numn);
%% Main loop for the analysis of each patch
[numm numn]
for cntm = 1:numm
    for cntn = 1:numn
        [cntm cntn]
        phi = phiwhole((1:paraSST.N)+(cntm-1)*paraSST.N/2,(1:paraSST.N)+(cntn-1)*paraSST.N/2);
        % Zero padding for non-periodic image
        [m n] = size(phi);
        dm = 0;
        dn = 0;
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
        
        % Search for peaks
        num_wave = 1;%determined by the hexegonal reference grain
        R_low = min(paraSST.R_low2,paraSST.R_low3)*(1-paraSST.spExtention); R_high = 2*max(paraSST.R_high2,paraSST.R_high3)*(1+paraSST.spExtention);
        R_low
        R_high
        %ss transform
        NB = [20 9];
        SPg = [20 20];
        ss_energy = SS_ct_polar(paraSST.num_direction,phi,SPg,NB,paraSST.rad,paraSST.is_real,R_low,R_high,paraSST.epsl,paraSST.red,paraSST.t_sc,paraSST.s_sc);
        valMax = max(ss_energy(:));
        ss_energy = ss_energy/valMax;
        szSST = size(ss_energy);
        pos1 = gdSzType:gdSzType:szSST(3)-gdSzType; pos2 = gdSzType:gdSzType:szSST(4)-gdSzType;
        ssEnergySpl = ss_energy(:,:,pos1,pos2);
        ssEnergySpl = permute(ssEnergySpl,[3,4,1,2]);
        ssEnergySpl = reshape(ssEnergySpl,[numel(pos1)*numel(pos2),szSST(1),szSST(2)]);
        adjMat = AdjMatRotInv(ssEnergySpl,sigma);
        [group_resid , svals, Lap , num_group_resid] = SpectralClustering_est(adjMat,3);
      %  group_resid
        typeMatCell{cntm,cntn} = zeros(numel(pos1)*numel(pos2),1);
        rec = [];
        tp = 0;
        for cnt = 1:num_group_resid
            posg = find(group_resid==cnt);
            if numel(posg) > sqrt(numel(pos1)*numel(pos2))
                tp = tp + 1;
                %posg
                rec = [rec,posg(1)];
                typeMatCell{cntm,cntn}(posg) = tp;
            end
        end
        typeMatCell{cntm,cntn} = reshape(typeMatCell{cntm,cntn},[numel(pos1),numel(pos2)]);
        typeCell{cntm,cntn} = ssEnergySpl(rec,:,:);
        ss_energy = reshape(ss_energy,[szSST(1),szSST(2),szSST(3)*szSST(4)]);
    end
end