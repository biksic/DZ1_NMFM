function [opti_omega, spekt_rad, ro]=sor_konvergencija(A)

tol=1e-8;
d=size(A); 
n=d(1);
x0=rand(n,1); %generiramo random vektor stupac 

x0=x0/norm(x0);

D=diag(diag(A)); 
R=triu(A)-D;  %gornjetrokutasta bez dijagonale
L=tril(A)-D;  %donjetrokutasta bez dijagonale

omega=[0:0.01:2]; %omega iz [0,2] korak 0,01


for k=1:201  % od 0 do 2 ima toliko omega ako brojimo i 0 (poèetak)
  T=(D+omega(k)*L)\((1-omega(k))*D-omega(k)*R);

  %spektralni radijus T je maskimalna svojstevna vrijednost matrice T
  [x,k_p,r,flag]=metoda_potencija(T,x0,tol);
  if(flag==0)
    ro(k)=max(abs(eig(T))); 
  else
   ro(k)=abs((x'*T*x)/(x'*x)); 
  end
  
  if (k==1)
    romin=ro(k);
    opti_omega=omega(1);
  
  else 
    if (ro(k)<romin)
      romin=ro(k);
      opti_omega=omega(k);
      
    end
    
  end
  
end


spekt_rad=romin; %samo zbog naziva

end