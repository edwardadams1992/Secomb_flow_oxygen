function rtf=rtflsp(x1,x2,xacc,hext,cext)
maxit=30;

fl=func(x1,hext,cext);
fh=func(x2,hext,cext);
if fl*fh>0
    rtf=-1;
    return    
end

if fl<0
    xl=x1;
    xh=x2;
else
    xl=x1;
    xh=x1;
    swap=fl;
    fl=fh;
    fh=swap;
end

dx=xh-xl;

for j=1:maxit
    rtf=xl+dx*fl/(fl-fh);
    f=func(rtf,hext,cext);
    if f<0
        del=xl-rtf;
        xl=rtf;
        fl=f;
    else
        del=xh-rtf;
        xh=rtf;
        fh=f;
    end
    dx=xh-xl;
    if abs(del)<xacc || f==0
        rtf=rtf;
        return
    end
    
end


rtf=-1;
return
end