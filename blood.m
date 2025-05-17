function [p,pp]=blood(c,h)

alphab=3.1e-5;
p50=27;%41.5
cs=0.5;
fn=3.0;%3.5%2.7; %n in hill equation
plow=0.1*p50;
phigh=5*p50;
clowfac=cs*(1-1/(1+(plow/p50)^fn));
chighfac=cs*(1-1/(1+(phigh/p50)^fn));
pphighfac=cs*fn/p50*((phigh/p50)^(fn-1))/((1+(phigh/p50)^fn)^2);


if h<1.e-6
    p=c/alphab;
    pp=alphab;
    return
end


if c<0
    p=c/alphab;
    pp=alphab;
    return
    
end

clow=clowfac*h+alphab*plow;

if c<clow
    p=c*plow/clow;
    pp=clow/plow;
    return
end

chigh=chighfac*h+alphab*phigh;

if c<chigh
    if c/h/cs<1
        ph=(((c/h/cs)/(1-c/h/cs))^(1/fn))*p50;
        pl=0;
    else
        ph=phigh;
        pl=plow;
        
    end
    
    hext=h;
    cext=c;
    p=rtflsp(pl,ph,0.001,hext,cext);
    
    if p<0
        p=rtbis(pl,ph,0.001,hext,cext);
    end
    pp=cs*h*fn/p50*((p/p50)^(fn-1))/((1+((p/p50)^fn))^2)+alphab;
    return

end

pphigh=pphighfac*h+alphab;
p=phigh+(c-chigh)/pphigh;
pp=pphigh;
return
end