function [q,kappa,lambda,nodpress,jj]=flow1(nseg,varyviscosity,constvisc,diam,lseg,varytargetshear,nodpress,ista,iend,flow_direction,known_flow_direction,ktau,kappa,solvetype,nnod,length_weight,hmat,kmat,nodtyp,nodseg,knowntyp,precond,nodelambda,presstarget1,nodeinflow,matrixdim,nodnod,lambda,nitmax2,totallength,eps,known_flow_weight)

kpress=1;
facfp=pi*1333/128/0.01*60/(1e6);%unit converter with pi
shearconstant=32/pi*(1e4)/60; %unit converter with pi

for iseg=1:nseg
    if varyviscosity==1
        viscosity=viscocr(diam(iseg),hd(iseg));
    else
        viscosity=constvisc;       
    end
    cond(iseg)=facfp*(diam(iseg)^4)/lseg(iseg)/viscosity; %C=1/R=Q/P=(pi/128)*D^4/(mu*L) - M
    shearfac(iseg)=shearconstant*viscosity/(diam(iseg)^3); %c_j
    
    %calc target shear
    if varytargetshear==1
        vesspress=((nodpress(ista(iseg)))+nodpress(iend(iseg)))/2; %vesspress=midpoint segment pressure=p_j
        vesspress=max([vesspress,10]);%vesspress=midpoint segment pressure=p_j
        vess_shear_target=100-86*exp(-5000*(log10(log10(vesspress)))^5.4); %(eq 13)
        sheartarget(iseg)=flow_direction(iseg)*vess_shear_target;
        
    else
        sheartarget(iseg)=flow_direction(iseg)*sheartarget1;
        
    end
    
    if known_flow_direction(iseg)==0
        hfactor1(iseg)=ktau*lseg(iseg)*shearfac(iseg)*kappa*cond(iseg)*sheartarget(iseg); %eq 11 - rhs
        hfactor2(iseg)=ktau*lseg(iseg)*((shearfac(iseg)*cond(iseg))^2);%eq12*ktau (ktau*Hik)
    else
        hfactor1(iseg)=known_flow_weight*shearfac(iseg)*kappa*cond(iseg)*sheartarget(iseg);
        hfactor2(iseg)=known_flow_weight*((shearfac(iseg)*cond(iseg))^2);
    end
end
    

if solvetype==2
        
[bvector,precond,kmat,hmat]=bvectorc(nnod,kpress,length_weight,hmat,kmat,nodtyp,nodseg,hfactor1,hfactor2,cond,knowntyp,ktau,precond,nodelambda,presstarget1,nodeinflow,nodpress,matrixdim,ista,iend,nodnod);      
for inod=1:matrixdim
    if inod<=nnod
       xvector(inod)=nodpress(inod);
    else
       xvector(inod)=lambda(inod-nnod);
    end
   
end

[xvector,jj]=sparsecgs(xvector,bvector,matrixdim,nitmax2,nnod,nodelambda,knowntyp,nodtyp,precond,kmat,hmat,nodnod,eps);

for inod=1:matrixdim

if inod<=nnod

nodpress(inod)=xvector(inod)*precond(inod);
else
lambda(inod-nnod)=xvector(inod)*precond(inod);
 
end


end

end
    
kappasum1=0;
kappasum2=0;

for iseg=1:nseg

q(iseg)=(nodpress(ista(iseg))-nodpress(iend(iseg)))*cond(iseg);

tau(iseg)=(nodpress(ista(iseg))-nodpress(iend(iseg)))*1333*diam(iseg)/lseg(iseg)/4;

kappasum2=kappasum2+lseg(iseg)*(tau(iseg)^2);
kappasum1=kappasum1+lseg(iseg)*abs(tau(iseg));


end
    
    
kappasum1=kappasum1/totallength;
kappasum2=kappasum2/totallength;
kappa=kappasum2/((kappasum1)^2);


%% chcekc conservation of flow

notconservedcount=0;

for inod=1:nnod
if knowntyp(inod)==1
currentflow=0;

    for i=1:nodtyp(inod)
    iseg=nodseg(i,inod);

    if ista(iseg)==inod

    currentflow=currentflow+q(iseg);
    else
    currentflow=currentflow-q(iseg);
    end
    end

    if abs(currentflow)>0.001
        notconservedcount=notconservedcount+1;
    end


end
end

%%
press_rms=0;
shear_rms=0;

for inod=1:nnod
press_rms=press_rms+length_weight(inod)*((nodpress(inod)-presstarget1)^2);
end

for iseg=1:nseg
shear_rms=shear_rms+lseg(iseg)*((tau(iseg)-sheartarget(iseg))^2);
end

press_rms=press_rms/totallength;
shear_rms=shear_rms/totallength;
total_dev=kpress*press_rms+ktau*shear_rms;
press_rms=sqrt(press_rms);
shear_rms=sqrt(shear_rms);
end