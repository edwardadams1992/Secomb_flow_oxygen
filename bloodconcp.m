function cvp=bloodconcp(p,h)

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
    cvp=alphab;
elseif p<plow
    cvp=clowfac*h/plow+alphab;
elseif p<phigh
    cvp=cs*h*fn/p50*((p/p50)^(fn-1))/((1+((p/p50)^fn))^2)+alphab;
else
    cvp=pphighfac*h+alphab;
end

end