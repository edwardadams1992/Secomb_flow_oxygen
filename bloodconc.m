function cv=bloodconc(p,h)

alphab=3.1e-5;
p50=27;%41.5
cs=0.5;
fn=3.0;%3.5%2.7; %n in hill equation
plow=0.1*p50;
phigh=5*p50;
clowfac=cs*(1-1/(1+(plow/p50)^fn));
chighfac=cs*(1-1/(1+(phigh/p50)^fn));
pphighfac=cs*fn/p50*((phigh/p50)^(fn-1))/((1+(phigh/p50)^fn)^2);


if p<0
    cv=alphab*p;
elseif p<plow
    cv=clowfac*h*p/plow+alphab*p;
elseif p<phigh
    cv=cs*h*(1-1/(1+(p/p50)^fn))+alphab*p; %eq 3/Q
else
    cv=(chighfac+(p-phigh)*pphighfac)*h+alphab*p;
end

end