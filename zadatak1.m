%load matrice1 
load matrice1.mat
format long;

%uvjetovanost matrica
uvjet_1=cond(A1);
uvjet_2=cond(A2);

%optimalan omega i spektralni radijus matrica
[op_om1,sp_rad1,ro_1]=sor_konvergencija(A1);
[op_om2,sp_rad2,ro_2]=sor_konvergencija(A2);

%spektralni radijus za A1
figure(1);
plot(omega,ro_1,'r-',op_om1,sp_rad1,'bx');
axis([0 2 0 1]);
xlabel('omega'); ylabel('spektralni radijus');
title('Graf spektralnih radijusa za A1');
legend('graf spektralnih radijusa','optimalan omega');
grid on;
print spektA1.pdf 

%spektraln radijus za A2
figure(2);
plot(omega,ro_2,'r-',op_om2,sp_rad2,'bx');
axis([0 2 0 1]);
xlabel('omega'); ylabel('spektralni radijus');
title('Graf spektralnih radijusa za A2');
legend('graf spektralnih radijusa','optimalan omega');
grid on;
print spektA2.pdf

%rješavanje sustava
[x_1,k_1,re_1,gr_1]=sor(A1,b1,x0,tol,op_om1);
[x_2,k_2,re_2,gr_2]=sor(A2,b2,x0,tol,op_om2);

%norme gresaka
figure(3);
semilogy(0:k_1,gr_1,'b-',0:k_2,gr_2, 'r-');
xlabel('broj iteracija k'); ylabel('greska ||x-x_k||');
legend('A1','A2');
title('Grafovi normi gresaka');
grid on;
print normegresaka1.pdf

%norme reziduala
figure(4);
semilogy(0:k_1,re_1,'b-',0:k_2, re_2,'r-');
xlabel('broj iteracija k'); ylabel('rezidual || r_k ||/|| b||');
legend('A1','A2');
title('Grafovi relativnih normi reziduala');
grid on;
print normereziduala1.pdf



