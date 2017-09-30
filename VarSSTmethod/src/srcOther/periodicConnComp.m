function [CC,colComponents] = periodicConnComp( binaryImg,is_periodic)
% binaryImg: 0 for the support of targets and 1 for space between targets
if nargin<2, is_periodic = 1; end;

CC = bwconncomp(1-binaryImg);
colComponents = zeros(size(binaryImg));
for i=1:CC.NumObjects
    colComponents(CC.PixelIdxList{i}) = i;
end
if is_periodic
    CC = ConnCompPeriodic(colComponents, CC);
    for i=1:CC.NumObjects
        colComponents(CC.PixelIdxList{i}) = i;
    end
end
end