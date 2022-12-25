function [x,k,re,flag]=metoda_potencija(A,x0,tol)
  
  k=0;
  x=x0/norm(x0); 
  re(1)=norm(A*x-(x'*A*x)*x);
  flag=1;
  
  while (re(k+1)>tol && flag==1) % flag dodan zbog napomenu 
    y=A*x;
    x=y/norm(y);
    k=k+1;
    re(k+1)=norm(A*x-(x'*A*x)*x); %racunamo reziduale
    if (k>100) %uvjet napomene
        flag=0; 
    end  
  end 

end