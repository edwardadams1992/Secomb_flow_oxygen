

function [viscosity]=viscor(d,h)

mcvcorr=1.187064456000345;
optw=1.1;
vplas=1.0466;
cpar=[0.8,-0.075,-11,12];
viscpar=[6,-0.085,3.2,-2.44,-0.06,0.645];
hdref=0.45;

dcorr=d*mcvcorr;

c=(cpar(1)+exp(cpar(2)*dcorr))*(-1+1/(1+10^(cpar(3))*(dcorr^cpar(4))))+1/(1+(10^cpar(3))*(dcorr^cpar(4)));

eta45=viscpar(1)*exp(viscpar(2)*dcorr)+viscpar(3)+viscpar(4)*exp(viscpar(5)*(dcorr^viscpar(6)));

hdfac=(((1-h)^c)-1)/(((1-hdref)^c)-1);

etarel=(1+(eta45-1)*hdfac*(dcorr/(dcorr-optw))^2)*(dcorr/(dcorr-optw))^2;

viscosity=etarel*vplas;

end