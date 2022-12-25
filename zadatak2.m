%load matrice2
load matrice2.mat
format long;

%raèunanje uvjetovanosti matrica
uvjet_3=cond(A3);
uvjet_4=cond(A4);

% hessenbergova forma matrica
[H_3,Q_3]=hessenberg(A3);
[H_4,Q_4]=hessenberg(A4);

%schur_qr
[T_3,U_3,ind_3,k3_]=schur_qr(H_3,tol);
[T_4,U_4,ind_4,k4_]=schur_qr(H_4,tol);


%gledamo samo dijagonalu jer smo matricu nastimali sami pa nema kompleksnih svojstvenih vrijednosti
sv_3=diag(T_3);
sv_4=diag(T_4);

%rješavanje sustava
[x_3,k_3,re_3,gr_3]=cg(A3,b3,x0,tol);
[x_4,k_4,re_4,gr_4]=cg(A4,b4,x0,tol);

%grafovi normi grešaka
figure(1);
semilogy(0:k_3,gr_3,'b-',0:k_4,gr_4, 'r-');
xlabel('broj iteracija k'); ylabel('greska ||x-x_k||');
legend('A3','A4');
title('Grafovi normi gresaka');
grid on;
print normegreska2.pdf

%grafovi relativnih normi reziduala
figure(2);
semilogy(0:k_3,re_3,'b-',0:k_4, re_4,'r-');
xlabel('broj iteracija k'); ylabel('|| r_k ||/|| b||');
legend('A3','A4');
title('Grafovi relativnih normi reziduala');
grid on;
print normereziduala2.pdf


