function hist = skeletonHist(skeleton)
sz = size(skeleton);
for cnt1 = 1:sz(3)
    for cnt2 = 1:sz(4)
        data = skeleton(:,:,cnt1,cnt2);
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
    end
end

