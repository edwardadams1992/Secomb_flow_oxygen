
function [nodpress,nodtyp,maxerr]=solve(nodsegm,nnod,nodtyp,nodseg,nnodbc,bcprfl,cond,bctyp,nodpress,omega,tol,nodnod,bcnod,nitmax)

wk=zeros(nodsegm,nnod);

for inod=1:nnod    
    if nodtyp(inod)==1
        wk(1,inod)=0;
    end
    
    if nodtyp(inod)>1
        condsum=0;
        for i=1:nodtyp(inod)
            iseg=abs(nodseg(i,inod));
            condsum=condsum+cond(iseg);
            wk(i,inod)=cond(iseg);
        end
        for i=1:nodtyp(inod)
            wk(i,inod)=wk(i,inod)/condsum;
        end
        
        
    end
    
end

for inodbc=1:nnodbc
    inod=bcnod(inodbc);
    if bctyp(inodbc)==0
        nodpress(inod)=bcprfl(inodbc);
        nodtyp(inod)=-1;
    else
        wk(1,inod)=bcprfl(inodbc)/cond(nodseg(1,inod));
    end
end


for niter=1:nitmax
    maxerr=0;
    for inod=1:nnod
        if nodtyp(inod)==1
            press1=omega*(nodpress(nodnod(1,inod))+wk(1,inod)-nodpress(inod));
        end
        if nodtyp(inod)>=2
            pcondsum=0;
            for i=1:nodtyp(inod)
                pcondsum=pcondsum+wk(i,inod)*nodpress(nodnod(i,inod));
                press1=omega*(pcondsum-nodpress(inod));
            end
        end
        
        if nodtyp(inod)>=1
            nodpress(inod)=nodpress(inod)+press1;
            if abs(press1)>=maxerr
                maxerr=abs(press1);
                
                errnode=inod;

                
            end
        end
              
    end
    
    if maxerr<tol
        
%         disp('linear converged')
        break;
    else
%         disp('linear not converged')
        
    end
end

for inod=1:nnod
    if nodtyp(inod)==-1
        nodtyp(inod)=1;
    end
end
end