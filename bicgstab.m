
function out=bicgstab(a,b,x,n,bigstaberr,itmax)
lu=0;
for i=1:n
    
    r(i)=0;
    for j=1:n
        r(i)=r(i)+a(i,j)*x(j);
    end
    r(i)=r(i)-b(i);
    p(i)=r(i);
    rs(i)=1;
    lu=lu+r(i)*rs(i);
    
end




kk=1;

while true
    t1=0;
    
    for i=1:n
        v(i)=0;
        for j=1:n
            v(i)=v(i)+a(i,j)*p(j);
        end
        t1=t1+v(i)*rs(i);
    end
    delta=-lu/t1;
    for i=1:n
        s(i)=r(i)+delta*v(i);
    end
    for i=1:n
        t(i)=0;
        for j=1:n
            t(i)=t(i)+a(i,j)*s(j);
        end
    end
    t1=0;
    t2=0;
    for i=1:n
        
        t1=t1+s(i)*t(i);
        t2=t2+t(i)*t(i);
    end
    gamma1=-t1/t2;
    err=0;
    lunew=0;
    
    for i=1:n
        r(i)=s(i)+gamma1*t(i);
        er=delta*p(i)+gamma1*s(i);
        x(i)=x(i)+er;
        if abs(er)>err
            err=abs(er);
        end
        lunew=lunew+r(i)*rs(i);
    end
    beta=lunew*delta/(lu*gamma1);
    lu=lunew;
    for i=1:n
        p(i)=r(i)+beta*(p(i)+gamma1*v(i));
    end
   
    kk=kk+1;

    if kk>itmax || err<bigstaberr
        out=x;
        disp(append('bigstab iterations =',num2str(kk)))
    disp(append('bigstab error =',num2str(err)))
        
        break
    end
    


    
end