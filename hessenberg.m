function [H,Q]=hessenberg(A)
  format long;
  n=length(A);
  Q=eye(n);

  if(n<=2) %vec u zadanoj normi
    Q=eye(n);
    H=A;
  else
    for k=1:n-2
      x=A(k+1:n,k);
      [v,beta]=gallery('house',x); %za Householderov reflektor
      
      A1=A(k+1:n,1:n);
      A(k+1:n,1:n)=A1-beta*v*(v'*A1);
     
      A2=A(1:n,k+1:n);
      A(1:n,k+1:n)=A2-beta*(A2*v)*v';
      
      A(k+2:n,k)=zeros(n-k-1,1); %poništili elemente stavljajuæi ih na nulu 
      
      Q(1:n,k+1:n) = Q(1:n,k+1:n)-beta*(Q(1:n,k+1:n)*v)*v';
      
    end
 H=A;   
end