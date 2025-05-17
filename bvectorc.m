

function [output,precond,kmat,hmat]=bvectorc(nnod,kpress,length_weight,hmat,kmat,nodtyp,nodseg,hfactor1,hfactor2,cond,knowntyp,ktau,precond,nodelambda,presstarget1,nodeinflow,nodpress,matrixdim,ista,iend,nodnod)

%rhs creation - p and lambda
for inod=1:nnod
    
    hmat(0+1,inod)=kpress*length_weight(inod);
    kmat(0+1,inod)=0;
    
    for i=1:nodtyp(inod)
        iseg=nodseg(i,inod);
        hmat(0+1,inod)=hmat(0+1,inod)+hfactor2(iseg); %eq 11 coeffcients. (k_p*w_i)+(sum over all connected segs to inod: ktau*H_ik)
        hmat(i+1,inod)=-hfactor2(iseg); %-H_ik at the ith connected seg
        kmat(0+1,inod)=kmat(0+1,inod)+cond(iseg)*max([1,1000*ktau]); % sum M over all connected segs to inod
        kmat(i+1,inod)=-cond(iseg)*max([1,1000*ktau]); % -M at the ith connected seg

        
    end
    
end

for inod=1:nnod %precond required for solving linear system
    if knowntyp(inod)~=0
        precond(inod)=1/sqrt(hmat(0+1,inod));
        
        if knowntyp(inod)~=3
            precond(nodelambda(inod))=sqrt(hmat(0+1,inod))/kmat(0+1,inod);
        end
    else
        precond(inod)=1;
        
    end
    
end

for inod=1:nnod
    if knowntyp(inod)==0 %if the pressure bc is known, give it the assigned pressure!
        output(inod)=nodpress(inod);
    else
        output(inod)=kpress*length_weight(inod)*presstarget1; %else give it length weighted pressure target p_0k
        
        if knowntyp(inod)~=3
            output(nodelambda(inod))=nodeinflow(inod)*max([1,1000*ktau]); %if not unknown bc - rhs lambda=nodeinflow
        end
        
        for i=1:nodtyp(inod)
            iseg=nodseg(i,inod);
            
            if ista(iseg)==inod
                output(inod)=output(inod)+hfactor1(iseg); %rhs(node)=sum t_oj*c_j*M_ji*lj over all connected inflow segs (rhs eq11)
            end
            if iend(iseg)==inod
                output(inod)=output(inod)-hfactor1(iseg); %rhs(node)=sum -(t_oj*c_j*M_ji*lj) over all connected outflow segs (rhs eq11)
            end
            currentnode=nodnod(i,inod);
            if knowntyp(currentnode)==0
                output(inod)=output(inod)-hmat(i+1,inod)*nodpress(currentnode);
                if knowntyp(inod)~=3
                    output(nodelambda(inod))=output(nodelambda(inod))-kmat(i+1,inod)*nodpress(currentnode);%if not unknown bc, calc lambda rhs: Qin(bc)-sum Q (M*p)
                end
            end
        end
        
        
    end
    
end


for i=1:matrixdim
    
    output(i)=output(i)*precond(i); %multiple by precondition for solving
end

end