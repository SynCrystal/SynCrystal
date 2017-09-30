function mat = skeletonPick(skeleton,refSkt,para)
% This code decide the type of crystal configuration according to the
% reference phase space skeleton.
%
% Input:
% skeleton: the phase space skeleton of the crystal image
% refSkt: the reference phase space skeleton
% para: a structure containing all parameters
method = para.typePick;
switch method
    case 1 % Eucliean distance
        numType = numel(refSkt);
        szSkt = size(skeleton);
        mat = zeros(szSkt(3),szSkt(4));
        for cnt = 1:numType
            temp = refSkt{cnt}(1:ceil(end/2),:);
            mxVal = sum(temp(:));
            if mxVal>para.clsThre/2
                temp = temp/mxVal;
            end
            refSkt{cnt} = temp;
         %   figure;imagesc(refSkt{cnt});axis image;colorbar;%caxis([0,1]);
        end
        %figure;
        for cnt2 = 1:szSkt(4)
            for cnt1 = 1:szSkt(3)
                temp = skeleton(1:ceil(end/2),:,cnt1,cnt2);
                mxVal = sum(temp(:));
                if mxVal>para.clsThre/2
                    temp = temp/mxVal;
                end
                % imagesc(temp);axis image;colorbar;caxis([0,1]);
                val = zeros(1,numType);
                for cnt = 1:numType
                    val(cnt) = norm(temp-refSkt{cnt});%sum(abs(temp(supp{cnt})-refVal{cnt}));
                end
                [~,pos] = min(val);
                mat(cnt1,cnt2) = pos;
                %  val
                %   pos
                %  pause;
            end
        end
    case 2 % 1D approximate EMD
        numType = numel(refSkt);
        szSkt = size(skeleton);
        mat = zeros(szSkt(3),szSkt(4));
        refHist = cell(1,numType);
        for cnt = 1:numType
            temp = refSkt{cnt};
            mxVal = sum(temp(:));
            if mxVal>para.clsThre/2
                temp = temp/mxVal;
            else
                temp = ones(size(temp))/size(temp,2);
            end
            refSkt{cnt} = temp;
            [~,pos] = max(sum(temp,2));
            refHist{cnt} = temp(pos,:)';
           % figure;imagesc(refSkt{cnt});axis image;colorbar;%caxis([0,1]);
        end
        %figure;
        for cnt2 = 1:szSkt(4)
            for cnt1 = 1:szSkt(3)
                temp = skeleton(:,:,cnt1,cnt2);
                mxVal = sum(temp(:));
                if mxVal>para.clsThre/2
                    temp = temp/mxVal;
                else
                    temp = ones(size(temp))/size(temp,2);
                end
                [~,pos] = max(sum(temp,2));
                tempHist = temp(pos,:)';
                % imagesc(temp);axis image;colorbar;caxis([0,1]);
                val = zeros(1,numType);
                for cnt = 1:numType
                    val(cnt) =  distVec(tempHist,refHist{cnt});% EMD
                end
                [~,pos] = min(val);
                mat(cnt1,cnt2) = pos;
                %  val
                %   pos
                %  pause;
            end
        end
    case 3 % 2D EMD
        numType = numel(refSkt);
        szSkt = size(skeleton);
        mat = zeros(szSkt(3),szSkt(4));
        refHist = cell(1,numType);
        for cnt = 1:numType
            temp = refSkt{cnt};
            mxVal = sum(temp(:));
            if mxVal>para.clsThre/5
                temp = temp/mxVal;
            else
                temp = ones(size(temp));
            end
            refSkt{cnt} = temp;
            [~,pos] = max(sum(temp,2));
            if pos == 1
                temp = temp(1:3,:);
                temp = temp/sum(temp(:));
                refHist{cnt} = temp;
            else
                temp = temp(pos-1:pos+1,:);
                temp = temp/sum(temp(:));
                refHist{cnt} = temp;
            end
           % figure;imagesc(refSkt{cnt});axis image;colorbar;%caxis([0,1]);
        end
        %figure;
        for cnt2 = 1:szSkt(4)
            for cnt1 = 1:szSkt(3)
                temp = skeleton(:,:,cnt1,cnt2);
                mxVal = sum(temp(:));
                if mxVal>para.clsThre/5
                    temp = temp/mxVal;
                else
                    temp = ones(size(temp));
                end
                [~,pos] = max(sum(temp,2));
                tempHist = temp(pos-1:pos+1,:);
                tempHist = tempHist/sum(tempHist(:));
                % imagesc(temp);axis image;colorbar;caxis([0,1]);
                val = zeros(1,numType);
                for cnt = 1:numType
                    val(cnt) =  EMD(tempHist,refHist{cnt});% EMD
                end
                [~,pos] = min(val);
                mat(cnt1,cnt2) = pos;
                %  val
                %   pos
                %  pause;
            end
        end
end
