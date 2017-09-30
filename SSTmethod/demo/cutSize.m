function mat2 = cutSize(mat,dimension,ratio)
[m,n] = size(mat);
if (dimension == 1)
    mnew = round(m*ratio);
    mat2 = mat(1:mnew,:);
end
if (dimension==2)
    nnew = round(n*ratio);
    mat2 = mat(:,1:nnew);
end