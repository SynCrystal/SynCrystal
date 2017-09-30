function Mout = skeletonShow(Min,num)
Mout = zeros(num,num);
sz = size(Min);
Mtemp = zeros(sz(1),num);
for cnt = 1:num
    pos = 1+floor(sz(2)*cnt/num);
    if pos == sz(2) + 1
        pos = pos -1;
    end
    Mtemp(:,cnt) = Min(:,pos);
end
for cnt = 1:num
    pos = 1+floor(sz(1)*cnt/num);
    if pos == sz(1) + 1
        pos = pos -1;
    end
    Mout(cnt,:) = Mtemp(pos,:);
end