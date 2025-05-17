
function [xvec,jj]=sparsecgs(xvector,bvector,matrixdim,nitmax2,nnod,nodelambda,knowntyp,nodtyp,precond,kmat,hmat,nodnod,eps)




b=bvector;
x=xvector;
n=matrixdim;

itmax=nitmax2;


r=zeros(1,n);
v=zeros(1,n);
p=zeros(1,n);



%% Amultiply
[r]=Amulti(x,nnod,nodelambda,knowntyp,nodtyp,precond,kmat,hmat,nodnod);

%%

for i=1:n
    r(i)=b(i)-r(i);
    p(i)=r(i);
end

rsold=dot(r,r);
jj=0;


while 1

    [v]=Amulti(p,nnod,nodelambda,knowntyp,nodtyp,precond,kmat,hmat,nodnod);

    pv=dot(p,v);
    

    
    for ii=1:n %issue here at j=15;
        x(ii)=x(ii)+(rsold/pv*p(ii));
        r(ii)=r(ii)-(rsold/pv*v(ii));
    end
    
    
    
    
    rsnew=dot(r,r);

    for ii=1:n
        p(ii)=r(ii)+rsnew/rsold*p(ii); 
    end
    
    jj=jj+1;
    rsold=rsnew;
    

 
    if rsnew<n*eps*eps || jj>itmax 

    disp(['cgsymm= ',num2str(rsnew)])
        break;
        
    end

end
xvec=x;
disp(['jj= ',num2str(jj)])
% disp(['cgsymm= ',num2str(rsnew)])

end
