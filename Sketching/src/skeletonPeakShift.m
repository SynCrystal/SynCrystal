function [scaleN,shiftAgl] = skeletonPeakShift(SST,avgdx,avgdy)
sz = size(SST);
scaleN = zeros(sz(3),sz(4));
shiftAgl= zeros(sz(3),sz(4));
for cnt1 = 1:sz(3)
    for cnt2 = 1:sz(4)
        data = SST(:,:,cnt1,cnt2);
        val = max(data(:));
        if val>0
            data = data/val;
        end
        pos = find(data<0.3);
        data(pos) = 0;
        %SST(:,:,cnt1,cnt2) = data;
        temp = [0;sum(data,2);0];
        [pks,locs] = findpeaks(temp);
        numPks = numel(pks);
        if numPks == 1
            [M,I] = max(data(:));
            [I_row, I_col] = ind2sub(size(data),I);
            R = sqrt(avgdx(I_row,I_col,cnt1,cnt2).^2+avgdy(I_row,I_col,cnt1,cnt2).^2);
            scaleN(cnt1,cnt2) = R;
            shiftAgl(cnt1,cnt2) = acos(avgdx(I_row,I_col,cnt1,cnt2)/scaleN(cnt1,cnt2));
        else if numPks >1
                ed = round((locs(1)+locs(2))/2)-1;
                data(ed:end,:) = 0;
                [M,I] = max(data(:));
                [I_row, I_col] = ind2sub(size(data),I);
                R = sqrt(avgdx(I_row,I_col,cnt1,cnt2).^2+avgdy(I_row,I_col,cnt1,cnt2).^2);
                scaleN(cnt1,cnt2) = R;
                shiftAgl(cnt1,cnt2) = acos(avgdx(I_row,I_col,cnt1,cnt2)/scaleN(cnt1,cnt2));
            else
                R = 1e10;
                scaleN(cnt1,cnt2) = 0;
                shiftAgl(cnt1,cnt2) = 0;
            end
        end
        %         if R<20
        %             R
        %             avgdx(I_row,I_col,cnt1,cnt2)
        %             avgdy(I_row,I_col,cnt1,cnt2)
        %             I_row
        %             I_col
        %         pause;
        %         end
    end
end

