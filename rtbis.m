function rtb=rtbis(x1,x2,xacc,hext,cext)

maxit=30;


f=func(x1,hext,cext);
fmid=func(x2,hext,cext);

if f*fmid>=0
    rtb=0;
    return
end

if f<0
    rtb=x1;
    dx=x2-x1;    
else
    rtb=x2;
    dx=x1-x2;    
end

for j=1:maxit
    dx=dx*0.5;
    xmid=rtb+dx;
    fmid=func(xmid,hext,cext);
    if fmid<=0
        rtb=xmid;
    end
    if abs(dx)<xacc || fmid==0
        rtb=rtb;
        return
    end
    
    
end

rtb=0;
return

end