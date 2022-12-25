%odabrao matricu 13x13, ideja iz zadatka 6
x0=rand(13,1);
x0=x0/norm(x0);
x=ones(13,1);
tol=1e-8;
n=13;

%jednostruke svojstevene vrijednosti koje nisu bliske po modulu
spektar3=[0.01:0.05:0.63];
spektar3(1)=100; 
spektar3(13)=500;
D3=diag(spektar3);
X3=rand(13,13);
[Q3,R3]=qr(X3);
A3=Q3*D3*Q3';
b3=A3*x;

%višestruke svojstvene vrijednosti, relativno bliske po modulu 
spektar4=[10 10 10 10 10 10 10 10 10 10 10.5 10.5 11 ];
D4=diag(spektar4);
X4=rand(13,13);
[Q4,R4]=qr(X4);
A4=Q4*D4*Q4';
b4=A4*x;
