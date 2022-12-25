%odabrao matricu 13x13, ideja iz zadatatka 6 
x0=rand(13,1);
x0=x0/norm(x0);
x=ones(13,1);
tol=1e-8;
omega=[0:0.01:2];

%jedna svojstevena vrijednost se razlikuje od ostatka
spektar1=(0.1:0.01:0.22);
spektar1(6)=100;
D1=diag(spektar1);
X1=rand(13,13);
[Q1,R1]=qr(X1); 
A1=Q1*D1*Q1';
b1=A1*x;

%sve svojstvene vrijednosti bliske
spektar2=(0.1:0.01:0.22);
D2=diag(spektar2);
X2=rand(13,13);
[Q2,R2]=qr(X2);
A2=Q2*D2*Q2';
b2=A2*x;
