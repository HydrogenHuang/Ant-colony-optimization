function ACO
%

m=3;
n=4;
alpha=1;
beta=2;
rho=0.5;
oo=9999;
d = [oo 3 1 2;
     3 oo 5 4;
     1 5 oo 2;
     2 4 2 oo];
 %Initailize
 Route = 1-eye(n,n);
 Length = 0;
 for i=1:n
     