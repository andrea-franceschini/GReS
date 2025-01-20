function mat_new = expandMat(mat)
s1 = size(mat,1);
s2 = size(mat,2);
mat_new(2*(1:s1)'-1,2*(1:s2)'-1) = mat;
mat_new(2*(1:s1)',2*(1:s2)') = mat;
end

