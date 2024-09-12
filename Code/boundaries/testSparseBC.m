function [tr,tc,td] = testSparseBC(A,dofs,method)
% Testing speed of BC application to sparse matrices
n = size(A,1);
d_dofs = (dofs-1)*n+dofs; % get lin indices of diagonal
switch method
    case 1 % working with indices
        tic 
        A(dofs,:)=0;
        tr = toc;
        tic
        A(:,dofs)=0;
        tc = toc;
        tic
        A(d_dofs) = 1;
        td = toc;
    case 2 % row indices with linear index
        tic
        [r,c] = find(A(dofs,:));
        A(sub2ind(size(A), dofs(r), c)) = 0;
        tr = toc;
        tic
        A(:,dofs)=0;
        tc = toc;
        tic
        A = A - diag(diag(A)) + speye(n);
        td = toc;
    case 3
        tic
        [r,c,s] = find(A(dofs,:));
        A = A - sparse(dofs(r), c, s, size(A,1), size(A,2));
        tr = toc;
        tic
        A(:,dofs)=0;
        tc = toc;
        tic
        A = A - diag(diag(A)) + speye(n);
        td = toc;
end
end

