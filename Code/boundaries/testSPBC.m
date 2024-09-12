n = 25000;
nd = n/10;
A = sprand(n,n,0.02);
dof = sort(randi([1 n],nd,1));
[tr1,tc1,td1]=testSparseBC(A,dof,1);
[tr2,tc2,td2]=testSparseBC(A,dof,2);
[tr3,tc3,td3]=testSparseBC(A,dof,3);