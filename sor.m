function [x,k,re,gr]=sor(A,b,x0,tol,omega)
  k=0;
  n=length(x0);
  x_rj=ones(n,1);
  re(1)=norm(b-A*x0)/norm(b);
  gr(1)=norm(x_rj-x0);
  
  while (norm(b-A*x0)/norm(b))>tol %algoritam za sor metodu 
    k=k+1;
    for i=1:n
      x0(i)=(1-omega)*x0(i);
      pom=b(i);
      for j=1:i-1
        pom=pom-A(i,j)*x0(j);
      end
      for j=i+1:n
        pom=pom-A(i,j)*x0(j);
      end
      x0(i)=x0(i)+pom*omega/A(i,i);
    end
    re(k+1)=norm(b-A*x0)/norm(b);
    gr(k+1)=norm(x_rj-x0);
  end
  x=x0;
end
