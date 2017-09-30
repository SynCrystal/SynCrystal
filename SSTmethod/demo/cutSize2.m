function mat2 = cutSize2(mat,pct)
[m,n] = size(mat);
range1 = round(m*pct):m-round(m*pct);
range2 = round(n*pct):n-round(n*pct);
mat2 = mat(range1,range2);