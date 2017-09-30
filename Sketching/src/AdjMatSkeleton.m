function [mat,vecCell] = AdjMatSkeleton(dataIn,sigma,method,thre)
if nargin < 3
    method = 1;
end
sz = size(dataIn);
mat = zeros(sz(3),sz(3));
%figure;

switch method
    case 1 % Euclidean distance of the skeleton
        pos = find(dataIn<0.1); % simple dimension reduction by thresholding
        dataIn(pos) = 0;
        for cnt = 1:sz(3)
            val = max(max(dataIn(:,:,cnt)));
            if val~=0
                dataIn(:,:,cnt) = dataIn(:,:,cnt)/val;
            end
            pos = find(dataIn(:,:,cnt)<thre);
            tempd = dataIn(:,:,cnt);
            tempd(pos) = 0;
            dataIn(:,:,cnt) = tempd; 
        end
        for cnt1 = 1:sz(3)-1
            for cnt2 = cnt1+1:sz(3)
                v1 = dataIn(:,:,cnt1);
                v2 = dataIn(:,:,cnt2);
                mat(cnt1,cnt2) = gauss(norm(v1(:)-v2(:))/sqrt(length(v1(:))),sigma);
            end
        end
        vecCell = [];
        
    case 2 % Earth Mover's distance of the Euclidean
        for cnt = 1:sz(3)
            val = sum(sum(dataIn(:,:,cnt)));
            if val~=0
                dataIn(:,:,cnt) = dataIn(:,:,cnt)/val;
            end
        end
        for cnt1 = 1:sz(3)-1
            for cnt2 = cnt1+1:sz(3)
                v1 = dataIn(:,:,cnt1);
                v2 = dataIn(:,:,cnt2);
                mat(cnt1,cnt2) = gauss(skeletonEMD(v1,v2),sigma);
            end
        end
        vecCell = [];
        
    case 3 % selected features
        vecCell = cell(1,sz(3));
        tempCell = cell(1,sz(3));
        L = 0;
        
        pos = find(dataIn<0.1); % simple dimension reduction by thresholding
        dataIn(pos) = 0;
        for cnt = 1:sz(3)
            val = max(max(dataIn(:,:,cnt)));
            if val~=0
                dataIn(:,:,cnt) = dataIn(:,:,cnt)/val;
            end
            pos = find(dataIn(:,:,cnt)<thre);
            tempd = dataIn(:,:,cnt);
            tempd(pos) = 0;
            %subplot(1,3,1);
            %imagesc(tempd);axis image;colorbar;
            temp = sum(tempd,2);
            %subplot(1,3,2);plot(temp);
            locs = find(temp>0);
            numPks = numel(locs);
            L = max(numPks,L);
            tempCell{cnt} = zeros(1,numPks);
            
            for cntp = 1:numPks
                temp2 = [0,tempd(cntp,:),0];
                [pks2,locs2] = findpeaks(temp2);
                pos1 = find(locs2==2); pos2 = find(locs2==sz(2)+1);
                if (numel(pos1)*numel(pos2))
                    pks2(pos2) = [];
                end
                tempCell{cnt}(cntp) = numel(pks2);
            end
            
            %     subplot(1,3,3);plot(temp2);
            %     display('numPks:');
            %     numPks
            %     tempCell{cnt}
            %     pause;
        end
        for cnt = 1:sz(3)
            vecCell{cnt} = zeros(1,L);
            numPks = numel(tempCell{cnt});
            vecCell{cnt}(1:numPks) = tempCell{cnt};
        end
        for cnt1 = 1:sz(3)-1
            for cnt2 = cnt1+1:sz(3)
                mat(cnt1,cnt2) = gauss(norm(vecCell{cnt1}-vecCell{cnt2})/sqrt(L),sigma);
            end
        end
        
    case 4
        pos = find(dataIn<0.1);
        dataIn(pos) = 0;
        vecCell = cell(1,sz(3));
        tempCell = cell(1,sz(3));
        L = 0;
        %figure;
        for cnt = 1:sz(3)
            val = max(max(dataIn(:,:,cnt)));
            if val~=0
                dataIn(:,:,cnt) = dataIn(:,:,cnt)/val;
            end
            pos = find(dataIn(:,:,cnt)<thre);
            tempd = dataIn(:,:,cnt);
            tempd(pos) = 0;
            %subplot(1,3,1);
            %imagesc(tempd);axis image;colorbar;
            temp = sum(tempd,2);
            temp = [0;temp;0];
            %subplot(1,3,2);plot(temp);
            [pks,locs] = findpeaks(temp);
            numPks = numel(pks);
            L = max(numPks,L);
            tempCell{cnt} = zeros(1,numPks);
            if numPks == 1
                temp2 = sum(tempd);
                temp2 = [0,temp2,0];
                [pks2,locs2] = findpeaks(temp2);
                pos1 = find(locs2==2); pos2 = find(locs2==sz(2)+1);
                if (numel(pos1)*numel(pos2))
                    pks2(pos2) = [];
                end
                tempCell{cnt} = numel(pks2);
            else
                for cnt2 = 1:numPks
                    if cnt2 == 1
                        st = 1; ed = round((locs(1)+locs(2))/2)-1;
                        temp2 = sum(tempd(st:ed,:));
                        temp2 = [0,temp2,0];
                        [pks2,locs2] = findpeaks(temp2);
                        pos1 = find(locs2==2); pos2 = find(locs2==sz(2)+1);
                        if (numel(pos1)*numel(pos2))
                            pks2(pos2) = [];
                        end
                        tempCell{cnt}(1) = numel(pks2);
                    else if cnt2 == numPks
                            st = round((locs(cnt2-1)+locs(cnt2))/2)-1; ed = sz(1)-1;
                            temp2 = sum(tempd(st:ed,:));
                            temp2 = [0,temp2,0];
                            [pks2,locs2] = findpeaks(temp2);
                            pos1 = find(locs2==2); pos2 = find(locs2==sz(2)+1);
                            if (numel(pos1)*numel(pos2))
                                pks2(pos2) = [];
                            end
                            tempCell{cnt}(end) = numel(pks2);
                        else
                            st = round((locs(cnt2-1)+locs(cnt2))/2)-1;
                            ed = round((locs(cnt2)+locs(cnt2+1))/2)-1;
                            temp2 = sum(tempd(st:ed,:));
                            temp2 = [0,temp2,0];
                            [pks2,locs2] = findpeaks(temp2);
                            pos1 = find(locs2==2); pos2 = find(locs2==sz(2)+1);
                            if (numel(pos1)*numel(pos2))
                                pks2(pos2) = [];
                            end
                            tempCell{cnt}(cnt2) = numel(pks2);
                        end
                    end
                end
            end
            
            %     subplot(1,3,3);plot(temp2);
            %     display('numPks:');
            %     numPks
            %     tempCell{cnt}
            %     pause;
        end
        for cnt = 1:sz(3)
            vecCell{cnt} = zeros(1,L);
            numPks = numel(tempCell{cnt});
            vecCell{cnt}(1:numPks) = tempCell{cnt};
        end
        for cnt1 = 1:sz(3)-1
            for cnt2 = cnt1+1:sz(3)
                mat(cnt1,cnt2) = gauss(norm(vecCell{cnt1}-vecCell{cnt2})/sqrt(L),sigma);
            end
        end
end
mat = mat'+ mat + eye(sz(3));

