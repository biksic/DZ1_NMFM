function [T,U,ind_2,k]=schur_qr(H,tol) 

  k=0;
  n=length(H);
  ind_2=0;
  ind=zeros(1,n);
  I=eye(n);
  U=I;
  pom=1;
  
%QR faktorizacija 
while(pom==1 && n>2)                                            
    pom=0;
    j=1;
    ind_1=ind;
    ind_1(1)=1; %e1
    l=1;
    ind_2=ind;
    
    for i=1:n-1
      %je li neki ispoddijagonalni element 0, jel H strogo Hessenbergova?
      if (abs(H(i+1,i))<=(tol)*(abs(H(i,i))+abs(H(i+1,i+1)))) %uvjet zadatka
        l=l+3;
        ind_1(l)=i;
        H(i+1,1)=0; 
      end
      
      if (abs(H(i+1,i))>(tol)*(abs(H(i,i))+abs(H(i+1,i+1))))  
      %dali trebamo još iteracija? 
        ind_2(j)=i;
        if (j>1 && (ind_2(j)-ind_2(j-1))==1)  
            pom=1;                                                            
        end
        j=j+1;
      end
     
      if (i>1 && H(i+1,i)==0 && abs(ind_1(l)-ind_1(l-1))>1)
      %idemo u rekurziju ako imamo 0 ispod dijagonale
         H1=H((ind_1(l-1)+1):i,(ind_1(l-1)+1):i);                   
         [~,U1,~,k1]=schur_qr(H1,tol); 
         k=k1+k;
         U2=I;
         U2((ind_1(l-1)+1):i,(ind_1(l-1)+1):i)=U1;                   
         U=U*U2;
         H=U2'*H*U2;
      end
      
     
      
      if (i==n-1 && l>1) 
      %rekurzija za posljednju 0 ispod dijagonale
         H1=H(ind_1(l)+1:n,ind_1(l)+1:n);
         [~,U1,~,k1]=schur_qr(H1,tol);
         k=k1+k; 
         U2=I;
         U2((ind_1(l)+1):n,(ind_1(l)+1):n)=U1;                        
         U=U*U2;
         H=U2'*H*U2;
       end
       
      if (abs(H(i+1,i))<=(tol)*(abs(H(i,i))+abs(H(i+1,i+1))))
        H(i+1,i)=0; 
      end
    
    end
    

    %iteracija QR faktorizacije (givensove rotacije)
    if (pom==1)                                                       
      Q=eye(n);
      for i=1:n-1
        G=eye(n);
        vekt=[H(i,i);H(i+1,i)];
        [G_pom,y]=planerot(vekt);
        c=G_pom(1,1);
        s=G_pom(1,2);
        G(i:i+1,i:i+1)=[c, -s; s, c];
        H=G'*H;
        Q=Q*G;
      end
      
      H=H*Q;
      U=U*Q;
      k=k+1;
    end
    
end

if (max(ind_2)>0)
  ind_2=ind_2(ind_2>0);
else 
  ind_2=0;
end

T=H;
%ispod sporedne dijagoanle nam ispadaju mali brojevi pa ih možemo staviti na nula
%da budemo u skladu s zadatkom iz knjige
T=triu(T,-1);

end