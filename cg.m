function [x,k,re,gr]=cg(A,b,x0,tol)
  k=0;
  r=b-A*x0;
  d=r;

  x=x0;
  n=length(x0);
  x_rj=ones(n,1);
  
  %normu greške dodali zbog crtanja 
  gr(1)=norm(x_rj-x0);
  r1=r'*r;
  r2=r1;
  re(1)=norm(r)/norm(b);
  
  while re(k+1)>tol %algoritam za cg 
    alfa=(r1)/(d'*A*d);
    x=x+alfa*d;
    r=r-alfa*A*d;
    r2=r'*r;
    beta=r2/r1;
    d=r+beta*d;
    k=k+1;
    re(k+1)=norm(r)/norm(b);
    gr(k+1)=norm(x_rj-x);
    r1=r2;
    end
end
