function [mtiss,mptiss]=tissrate3(tissparam,nsp,p,prpoints)

for isp=1:nsp
    mtiss(isp)=0;
    mptiss(isp)=0;
end

if nsp==1
%oxygen
if prpoints==0
m0=tissparam(1,1);
else

        m0=4*tissparam(1,1);

end
pcr=tissparam(2,1);
po2=max([p(1),0]);
mtiss(1)=-m0*po2/(po2+pcr);
if po2>0
mptiss(1)=-m0*pcr/((po2+pcr)^2);
end
end


if nsp==2
%oxygen
if prpoints==0
m0=tissparam(1,1);
else
m0=4*tissparam(1,1);    
end
pcr=tissparam(2,1);
po2=max([p(1),0]);
mtiss(1)=-m0*po2/(po2+pcr);
if po2>0
mptiss(1)=-m0*pcr/((po2+pcr)^2);
end
%solute 1
if prpoints==0
m0=tissparam(1,2);
else
m0=5*tissparam(1,2);    
end
pcr=tissparam(2,2);
po2=max([p(2),0]);
mtiss(2)=-m0*po2/(po2+pcr);
if po2>0
mptiss(2)=-m0*pcr/((po2+pcr)^2);
end

end

end
